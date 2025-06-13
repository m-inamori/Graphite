#include <algorithm>
#include <cassert>
#include "../include/VCFIsolated.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFIsolated ////////////////////

VCFIsolated::VCFIsolated(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFRecord *>& rs, size_t nis,
							const Map& m, bool mg) :
									VCFBase(h, s), VCFImputable(m),
									records(rs), num_imputed_samples(nis),
									modify_genotypes(mg) { }

VCFIsolated::~VCFIsolated() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

vector<Haplotype> VCFIsolated::collect_haplotype_from_refs() const {
	vector<Haplotype>	haps;
	for(size_t i = num_imputed_samples; i < num_samples(); ++i) {
		haps.push_back(clip_haplotype(i, 0));
		haps.push_back(clip_haplotype(i, 1));
	}
	return haps;
}

vector<Haplotype> VCFIsolated::collect_haplotypes_mat(
										size_t sample_index) const {
	return collect_haplotype_from_refs();
}

vector<Haplotype> VCFIsolated::collect_haplotypes_pat(
										size_t sample_index) const {
	return collect_haplotype_from_refs();
}

VCFIsolated *VCFIsolated::divide_by_positions(size_t first, size_t last) const {
	vector<VCFRecord *>	sub_records(records.begin() + first,
									records.begin() + last);
	return new VCFIsolated(get_header(), get_samples(), sub_records,
							num_imputed_samples, get_map(), modify_genotypes);
}

vector<HaplotypePair> VCFIsolated::impute_cM(
									const vector<HaplotypePair>& prev_haps) {
	vector<HaplotypePair>	haps;
	for(size_t i = 0; i < prev_haps.size(); ++i) {
		const auto	prev_hap = prev_haps[i];
		const auto	hap = this->impute_cM_each_sample(prev_hap, i, true,
															modify_genotypes);
		haps.push_back(hap);
	}
	return haps;
}

void VCFIsolated::impute() {
	const Haplotype	h = Haplotype::default_value();
	vector<HaplotypePair>	haps(num_imputed_samples, make_pair(h, h));
	vector<VCFIsolated *>	vcf_cMs = VCFImputable::divide_by_cM(this);
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFIsolated	*vcf_cM = *p;
		haps = vcf_cM->impute_cM(haps);
		vcf_cM->records.clear();
		delete vcf_cM;
	}
	cout << num_imputed_samples << " samples are imputed." << endl;
}

VCFIsolated *VCFIsolated::create(const VCFSmall *orig_vcf,
											const VCFSmall *imputed_vcf,
											const STRVEC& samples,
											const STRVEC& references,
											const OptionSmall& op) {
	const auto	sample_columns = orig_vcf->extract_columns(samples);
	const auto	ref_columns = imputed_vcf->extract_columns(references);
	STRVEC	new_samples;
	for(auto p = sample_columns.begin(); p != sample_columns.end(); ++p) {
		new_samples.push_back(orig_vcf->get_samples()[*p-9]);
	}
	new_samples.insert(new_samples.end(), references.begin(), references.end());
	const auto	header = orig_vcf->trim_header(new_samples);
	vector<VCFRecord *>	records;
	// Create VCF first and add Record later to have samples in VCF
	VCFIsolated	*vcf = new VCFIsolated(header, new_samples, records,
											samples.size(), op.map, true);
	for(size_t i = 0; i < orig_vcf->size(); ++i) {
		VCFRecord	*record = orig_vcf->get_record(i);
		VCFRecord	*imputed_record = imputed_vcf->get_record(i);
		const STRVEC&	v0 = record->get_v();
		STRVEC	v(v0.begin(), v0.begin() + 9);
		for(auto q = sample_columns.begin(); q != sample_columns.end(); ++q)
			v.push_back(record->get_v()[*q]);
		for(auto q = ref_columns.begin(); q != ref_columns.end(); ++q)
			v.push_back(imputed_record->get_v()[*q]);
		VCFRecord	*new_record = new VCFRecord(v, vcf->get_samples());
		vcf->add_record(new_record);
	}
	return vcf;
}
