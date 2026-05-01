#include "../include/VCFGeno.h"
#include "../include/LargeSelfFamilyRef.h"
#include "../include/VCFSelfFillable.h"
#include "../include/VCFSelfImputable.h"
#include "../include/VCFSelfParentImputed.h"
#include "../include/VCFSelfProgenyImputed.h"
#include "../include/SampleManager.h"
#include "../include/KnownFamily.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// LargeSelfFamilyRef ////////////////////

VCFSelfParentImputed *LargeSelfFamilyRef::create_parent_phased_vcf(
													const KnownFamily *family,
													const VCFSmall *orig_vcf,
													const VCFGeno *phased_vcf,
													const Map& gmap,
													const Option& op) {
	const auto&	fsamples = family->get_samples();
	const auto	ids = phased_vcf->extract_columns(fsamples);
	vector<string>	samples(1, family->get_mat());
	for(size_t i = 1; i < fsamples.size(); ++i) {
		if(ids[i] != string::npos)
			samples.push_back(fsamples[i]);
	}
	if(samples.size() == 1)
		return NULL;
	
	vector<string>	prog_samples(samples.begin() + 1, samples.end());
	auto	vcf_progs = VCFGeno::extract_samples(prog_samples, orig_vcf);
	auto	records = RefCommon::merge_records(phased_vcf, vcf_progs, samples);
	return new VCFSelfParentImputed(samples, records, gmap, 0.01, orig_vcf);
}

VCFSelfProgenyImputed *LargeSelfFamilyRef::create_progeny_phased_vcf(
													const KnownFamily *family,
													const VCFSmall *orig_vcf,
													const VCFGeno *ref_vcf,
													const VCFGeno *phased_vcf,
													const Map& gmap,
													const Option& op) {
	const STRVEC&	samples = family->get_samples();
	const size_t	N = samples.size();
	const auto	ids = phased_vcf->extract_columns(samples);
	vector<size_t>	prog_ids;
	for(size_t i = 1; i < N; ++i) {
		if(ids[i] != string::npos)
			prog_ids.push_back(i);
	}
	
	const size_t	phased_prog_id = prog_ids[0];
	const string&	phased_sample = samples[phased_prog_id];
	// # Keep only the first progeny
	const set<size_t>	set_prog_ids(prog_ids.begin() + 1, prog_ids.end());
	STRVEC	new_samples;
	for(size_t i = 0; i < N; ++i) {
		if(set_prog_ids.find(i) == set_prog_ids.end())
			new_samples.push_back(samples[i]);
	}
	
	STRVEC	non_phased_samples;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p != phased_sample)
			non_phased_samples.push_back(*p);
	}
	
	const auto	*non_phased_vcf = VCFGeno::extract_samples(non_phased_samples,
																	orig_vcf);
	const auto	records = RefCommon::merge_records(phased_vcf, non_phased_vcf,
																new_samples);
	const auto	ref_haps = ref_vcf->create_ref_haps();
	return new VCFSelfProgenyImputed(new_samples, records, ref_haps,
											prog_ids[0], gmap, 0.01, orig_vcf);
}

VCFGeno *LargeSelfFamilyRef::extract_parents(
									const vector<VCFSelfFillable *>& vcfs) {
	STRVEC	samples;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		samples.push_back((*p)->get_samples()[0]);
	}
	
	const size_t	M = vcfs[0]->size();
	vector<GenoRecord *>	records;
	for(size_t i = 0; i < M; ++i) {
		const ll	pos = vcfs[0]->get_record(i)->get_pos();
		vector<int>	geno;
		for(size_t j = 0; j < vcfs.size(); ++j) {
			geno.push_back(vcfs[j]->get_record(i)->get_geno(0));
		}
		auto	record = new GenoRecord(pos, geno);
		records.push_back(record);
	}
	return new VCFGeno(samples, records, vcfs[0]->get_ref_vcf());
}

bool LargeSelfFamilyRef::is_intersect(const vector<string>& v,
											const set<string>& s) {
	for(auto p = v.begin(); p != v.end(); ++p) {
		if(s.find(*p) != s.end())
			return true;
	}
	return false;
}

VCFGeno *LargeSelfFamilyRef::impute(
								const vector<const KnownFamily *>& families,
								const VCFSmall *orig_vcf,
								VCFGeno *merged_vcf,
								const VCFGeno *ref_vcf,
								const Map& gmap, const Option& op) {
	if(families.empty())
		return NULL;
	
	const auto	*phased_vcf = VCFGeno::join(merged_vcf, ref_vcf,
												orig_vcf->get_samples());
	
	vector<VCFSelfImputable *>	vcfs;
	const STRVEC&	phased_samples = phased_vcf->get_samples();
	set<string>	set_phased(phased_samples.begin(), phased_samples.end());
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		STRVEC	samples = family->get_samples();
		if(set_phased.find(samples[0]) != set_phased.end()) {
			auto	*vcf1 = create_parent_phased_vcf(family, orig_vcf,
														phased_vcf, gmap, op);
			if(vcf1 != NULL)
				vcfs.push_back(vcf1);
		}
		else if(is_intersect(samples, set_phased)) {
			auto	*vcf2 = create_progeny_phased_vcf(family, orig_vcf, ref_vcf,
														phased_vcf, gmap, op);
			vcfs.push_back(vcf2);
		}
	}
	
	VCFSelfImputable::impute_VCFs(vcfs, op.num_threads);
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	if(merged_vcf != NULL)
		vcfs1.push_back(merged_vcf);
	
	if(vcfs1.empty())
		return NULL;
	
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
