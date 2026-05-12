#include <algorithm>
#include <cassert>
#include "../include/VCFGeno.h"
#include "../include/SelfFamily.h"
#include "../include/VCFSelfParentImputed.h"
#include "../include/VCFSelfProgenyImputed.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// SelfFamily ////////////////////

GenoRecord *SelfFamily::create_record(const GenoRecord *record,
										vector<size_t>& indices) {
	const int	pos = record->get_pos();
	vector<int>	geno(indices.size());
	for(size_t i = 0; i < indices.size(); ++i) {
		geno[i] = record->get_geno(indices[i]);
	}
	return new GenoRecord(pos, geno);
}

vector<size_t> SelfFamily::collect_progeny_indices(const Family *family,
										const set<string>& imputed_samples) {
	vector<size_t>	prog_indices;
	const auto&	progs = family->get_progenies();
	for(size_t j = 0; j < progs.size(); ++j) {
		const string&	prog = progs[j]->get_name();
		if(imputed_samples.find(prog) != imputed_samples.end())
			prog_indices.push_back(j);
	}
	return prog_indices;
}

VCFSelfImputable *SelfFamily::create_family_vcf(const VCFSmall *orig_vcf,
											const vector<GenoRecord *>& records,
											const vector<vector<int>>& ref_haps,
											const KnownFamily *family,
											const set<string>& imputed_samples,
											const OptionSmall& op) {
	// indices of imputed progenies
	const auto	prog_indices = collect_progeny_indices(family, imputed_samples);
	const auto&	mat = family->get_mat();
	const auto&	progs = family->get_progenies();
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	if(set_imputed_samples.find(mat) != set_imputed_samples.end()) {
		// If the parent is imputed, exclude imputed progenies from the VCF
		vector<size_t>	indices(1, 0);
		STRVEC	samples { mat };
		for(size_t i = 0; i < progs.size(); ++i) {
			const auto&	s = progs[i]->get_name();
			if(set_imputed_samples.find(s) == set_imputed_samples.end()) {
				samples.push_back(s);
				indices.push_back(i + 1);
			}
		}
		
		vector<GenoRecord *>	new_records;
		for(auto p = records.begin(); p != records.end(); ++p) {
			new_records.push_back(create_record(*p, indices));
		}
		
		return new VCFSelfParentImputed(samples, new_records,
												op.map, 0.01, orig_vcf);
	}
	else if(!prog_indices.empty()) {
		// If the parent is not imputed,
		// include only one imputed progeny in the VCF
		vector<string>	samples(1, mat);
		vector<size_t>	indices(1, 0);
		samples.push_back(progs[prog_indices[0]]->get_name());
		indices.push_back(prog_indices[0] + 1);
		for(size_t i = 0; i < progs.size(); ++i) {
			const auto&	s = progs[i]->get_name();
			if(set_imputed_samples.find(s) == set_imputed_samples.end()) {
				samples.push_back(s);
				indices.push_back(i + 1);
			}
		}
		
		vector<GenoRecord *>	new_records;
		for(auto p = records.begin(); p != records.end(); ++p) {
			new_records.push_back(create_record(*p, indices));
		}
		
		return new VCFSelfProgenyImputed(samples, new_records, ref_haps,
													0, op.map, 0.01, orig_vcf);
	}
	else {
		return NULL;
	}
}

VCFGeno *SelfFamily::impute(const VCFSmall *orig_vcf,
							const VCFGeno *imputed_vcf,
							const vector<vector<int>>& ref_haps,
							const vector<const KnownFamily *>& families,
							const vector<string>& imputed_samples,
							const OptionSmall& op) {
	const size_t	N = families.size();
	if(N == 0)
		return NULL;
	
	vector<VCFSelfImputable *>	vcfs;
	vector<VCFGenoBase *>	vcf_garbage;
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFGeno::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		auto	*vcf1 = create_family_vcf(orig_vcf, vcf->get_records(),
													ref_haps, family,
													set_imputed_samples, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
		}
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	VCFSelfImputable::impute_VCFs(vcfs, op.num_threads);
	
	cout << vcfs.size() << " self families have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
