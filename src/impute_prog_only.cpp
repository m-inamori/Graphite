#include <algorithm>

#include "../include/VCFHeteroHomo.h"
#include "../include/VCFHomoHomo.h"
#include "../include/VCFHeteroHeteroLite.h"
#include "../include/VCFJunkRecord.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/impute_prog_only.h"
#include "../include/ClassifyRecord.h"
#include "../include/option.h"

using namespace std;


//////////////////// ImputeProgOnly ////////////////////

void ImputeProgOnly::classify_record(size_t i, VCFFamilyRecord *record,
							const STRVEC& samples, const TypeDeterminer *td,
							std::vector<VCFHeteroHomoRecord *>& heho_records,
							std::vector<VCFImpFamilyRecord *>& other_records) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const auto	pair1 = CR->classify(record, td, false);
	const ParentComb	pc = pair1.first;
	const WrongType	wrong_type = pair1.second;
	const ll	pos = record->get_pos();
	const auto&	geno = record->get_geno();
	if(TypeDeterminer::is_homohomo(pc)) {
		auto	*record_ = new VCFHomoHomoRecord(pos, geno, i, wrong_type, pc);
		other_records[i] = record_->impute();
	}
	else if(TypeDeterminer::is_heterohomo(pc)) {
		heho_records[i] = new VCFHeteroHomoRecord(pos, geno, i, wrong_type, pc);
	}
	else if(pc == ParentComb::P01x01) {
		other_records[i] = new VCFHeteroHeteroLiteRecord(pos, geno, i,
																wrong_type, pc);
	}
	else {	// no candidate
		other_records[i] = new VCFJunkRecord(pos, geno, i, wrong_type);
	}
}

pair<vector<VCFHeteroHomoRecord *>, vector<VCFImpFamilyRecord *>>
ImputeProgOnly::classify_records(VCFFamily *vcf,
									const STRVEC& samples, Option *option) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const TypeDeterminer	*td = CR->get_TypeDeterminer(samples.size()-2,
																option->ratio);
	const size_t	M = vcf->size();
	vector<VCFHeteroHomoRecord *>	records1(M, NULL);
	vector<VCFImpFamilyRecord *>	records2(M, NULL);
	for(size_t i = 0; i < M; ++i) {
		VCFFamilyRecord	*record = vcf->get_family_record(i);
		classify_record(i, record, samples, td, records1, records2);
	}
	
	vector<VCFHeteroHomoRecord *>	heho_records;
	vector<VCFImpFamilyRecord *>	other_records;
	for(size_t i = 0; i < M; ++i) {
		auto	*r1 = records1[i];
		auto	*r2 = records2[i];
		if(r1 != NULL)
			heho_records.push_back(r1);
		else if(r2 != NULL)
			other_records.push_back(r2);
	}
	return make_pair(heho_records, other_records);
}

VCFRecord *ImputeProgOnly::fill_NA(const VCFRecord *record1,
											const STRVEC& samples) {
	const size_t	NA_len = samples.size() - record1->num_samples();
	const STRVEC&	v1 = record1->get_v();
	STRVEC	v = v1;
	v[8] = "GT";
	for(size_t c = 9; c < v.size(); ++c)
		v[c] = v[c].substr(0, 3);
	for(size_t i = 0; i < NA_len; ++i)
		v.push_back("./.");
	return new VCFRecord(v, samples);
}

VCFRecord *ImputeProgOnly::merge_record(const VCFRecord *record1,
										const VCFRecord *record2,
										const STRVEC& samples) {
	// TODO: convert VCFRecord to GenoRecord
	STRVEC	v = record1->get_v();
	const STRVEC&	v2 = record2->get_v();
	v[8] = "GT";
	for(auto p = v2.begin() + 9; p != v2.end(); ++p)
		v.push_back(p->substr(0, 3));
	return new VCFRecord(v, samples);
}

VCFSmall *ImputeProgOnly::merge_parents_progenies(const VCFSmall *vcf_parents,
												  const VCFSmall *vcf_progenies,
												  const STRVEC& samples) {
	const auto&	header = vcf_parents->trim_header(samples);
	vector<VCFRecord *>	records;
	const vector<VCFRecord *>&	parents_records = vcf_parents->get_records();
	size_t	j = 0;
	for(size_t i = 0; i < parents_records.size(); ++i) {
		const VCFRecord	*record1 = parents_records[i];
		if(j == vcf_progenies->size()) {
			records.push_back(fill_NA(record1, samples));
		}
		else {
			const VCFRecord	*record2 = vcf_progenies->get_record(j);
			if(record1->pos() == record2->pos()) {
				records.push_back(merge_record(record1, record2, samples));
				j += 1;
			}
			else {
				records.push_back(fill_NA(record1, samples));
			}
		}
	}
	return new VCFSmall(header, samples, records);
}

VCFHeteroHomoPP *ImputeProgOnly::merge_vcf(
							array<vector<VCFFillableRecord *>, 4>& rss,
							const STRVEC& samples, const Map& gmap,
							const VCFSmall *vcf) {
	vector<VCFFillableRecord *>	records;
	for(auto p = rss.begin(); p != rss.end(); ++p) {
		records.insert(records.end(), p->begin(), p->end());
	}
	std::sort(records.begin(), records.end(), 
				[](const VCFFillableRecord *lh, const VCFFillableRecord *rh)
				{ return lh->get_pos() < rh->get_pos(); });
	return new VCFHeteroHomoPP(samples, records, gmap, vcf);
}

VCFHeteroHomoPP *ImputeProgOnly::impute_prog_vcf_chr(
										const VCFSmall *parent_vcf,
										const VCFSmall *prog_vcf,
										const Map& gmap, const Option *option) {
	cout << "chr: " << parent_vcf->get_record(0)->chrom()
					<< parent_vcf->size() << " records "
					<< prog_vcf->size() << " records" << endl;
	// The assumption is that parent_vcf contains data only for the parents,
	// while prog_vcf contains data only for the progeny.
	auto	samples = parent_vcf->get_samples();
	const VCFSmall	*orig_vcf = merge_parents_progenies(parent_vcf,
														prog_vcf, samples);
	const auto&	prog_samples = prog_vcf->get_samples();
	samples.insert(samples.end(), prog_samples.begin(), prog_samples.end());
	auto	*merged_vcf = VCFHeteroHomoPP::merge(parent_vcf, prog_vcf,
											orig_vcf, samples, gmap, option);
	// Type conversion is necessary
	// before passing the data to the classify_records function.
	const auto&	rs = merged_vcf->get_records();
	vector<VCFFamilyRecord *>	records(rs.begin(), rs.end());
	auto	rss = VCFHeteroHomoPP::classify_records(samples, records, orig_vcf);
	const auto	header = parent_vcf->trim_header(samples);
	auto	*mat_vcf = new VCFHeteroHomoPP(samples,
										rss[static_cast<int>(FillType::MAT)],
										gmap, orig_vcf);
	auto	*pat_vcf = new VCFHeteroHomoPP(samples,
										rss[static_cast<int>(FillType::PAT)],
										gmap, orig_vcf);
	mat_vcf->impute();
	pat_vcf->impute();
	auto	*new_merged_vcf = merge_vcf(rss, samples, gmap, orig_vcf);
	new_merged_vcf->fill();
	return new_merged_vcf;
}
