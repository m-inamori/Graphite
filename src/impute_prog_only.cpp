#include <algorithm>

#include "../include/impute_prog_only.h"
#include "../include/VCFBothParentImputed.h"
#include "../include/option.h"

using namespace std;


//////////////////// ImputeProgOnly ////////////////////

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

VCFGenoBase *ImputeProgOnly::impute_prog_vcf_chr(
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
	auto	*merged_vcf = VCFGeno::convert(orig_vcf);
	
	vector<VCFFamilyRecord *>	records(merged_vcf->size(), NULL);
	for(size_t i = 0; i < merged_vcf->size(); ++i) {
		const GenoRecord	*record = merged_vcf->get_record(i);
		vector<int>	geno = record->get_geno();
		records[i] = new VCFFamilyRecord(record->get_pos(), geno);
	}
	auto	*vcf = new VCFBothParentImputed(merged_vcf->get_samples(),
												records, gmap, 0.01, orig_vcf);
	vcf->impute();
	return vcf->extract_by_samples(prog_samples);
}
