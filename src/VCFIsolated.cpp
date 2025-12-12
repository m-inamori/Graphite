#include <algorithm>
#include <cassert>
#include "../include/VCFIsolated.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFIsolated ////////////////////

VCFIsolated::VCFIsolated(const STRVEC& s, size_t nis, vector<GenoRecord *>& rs,
											const Map& m, const VCFSmall *vcf) :
									VCFClippable(s, m, vcf),
									records(rs), num_imputed_samples(nis) { }

VCFIsolated::~VCFIsolated() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

vector<VCFIsolated *> VCFIsolated::divide_by_cM() const {
	// If more than 1 cM but not more than 10 pieces and not more than 10 cM
	// are considered to be a block
	std::vector<VCFIsolated *>	vcfs;
	vector<GenoRecord *>	rs(1, records[0]);
	for(std::size_t i = 1; i < this->size(); ++i) {
		GenoRecord	*record = records[i];
		if(this->is_block(record, rs)) {
			rs.push_back(record);
		}
		else {
			vcfs.push_back(new VCFIsolated(get_samples(), num_imputed_samples,
												rs, get_map(), get_ref_vcf()));
			rs.clear();
			rs.push_back(record);
		}
	}
	vcfs.push_back(new VCFIsolated(get_samples(), num_imputed_samples,
												rs, get_map(), get_ref_vcf()));
	return vcfs;
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

vector<HaplotypePair> VCFIsolated::impute_cM(
									const vector<HaplotypePair>& prev_haps) {
	vector<HaplotypePair>	haps;
	for(size_t i = 0; i < prev_haps.size(); ++i) {
		const auto	prev_hap = prev_haps[i];
		const auto	hap = this->impute_cM_each_sample(prev_hap, i, true);
		haps.push_back(hap);
	}
	return haps;
}

void VCFIsolated::impute() {
	const Haplotype	h = Haplotype::default_value();
	vector<HaplotypePair>	haps(num_imputed_samples, make_pair(h, h));
	vector<VCFIsolated *>	vcf_cMs = this->divide_by_cM();
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFIsolated	*vcf_cM = *p;
		haps = vcf_cM->impute_cM(haps);
		vcf_cM->records.clear();
		delete vcf_cM;
	}
	cout << num_imputed_samples << " samples are imputed." << endl;
}

VCFIsolated *VCFIsolated::create(const VCFSmall *orig_vcf,
											const VCFGeno *imputed_vcf,
											const STRVEC& samples,
											const STRVEC& references,
											const OptionSmall& op) {
	const auto	columns = orig_vcf->extract_columns(samples);
	const auto	ref_columns = imputed_vcf->extract_columns(references);
	STRVEC	new_samples;
	for(auto p = columns.begin(); p != columns.end(); ++p) {
		new_samples.push_back(orig_vcf->get_samples()[*p-9]);
	}
	new_samples.insert(new_samples.end(), references.begin(), references.end());
	vector<GenoRecord *>	records;
	for(size_t i = 0; i < orig_vcf->size(); ++i) {
		VCFRecord	*record = orig_vcf->get_record(i);
		GenoRecord	*imputed_record = imputed_vcf->get_record(i);
		const ll	pos = imputed_record->get_pos();
		vector<int>	geno;
		for(auto p = columns.begin(); p != columns.end(); ++p)
			geno.push_back(Genotype::all_gt_to_int(record->get_v()[*p]));
		for(auto q = ref_columns.begin(); q != ref_columns.end(); ++q)
			geno.push_back(imputed_record->get_geno()[*q]);
		GenoRecord	*new_record = new GenoRecord(pos, geno);
		records.push_back(new_record);
	}
	return new VCFIsolated(new_samples, samples.size(),
									records, op.map, orig_vcf);
}
