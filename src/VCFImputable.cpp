#include "../include/VCFImputable.h"

using namespace std;

vector<int> VCFImputable::get_int_gts(size_t sample_index) const {
	vector<int>	int_gts;
	for(size_t i = 0; i < this->size(); ++i)
		int_gts.push_back(get_record(i)->unphased(sample_index));
	return int_gts;
}

Haplotype VCFImputable::clip_haplotype(size_t sample_index, int side) const {
	const auto	hap = this->clip_raw_haplotype(sample_index, side);
	return Haplotype(hap, sample_index, side);
}

HaplotypePair VCFImputable::impute_cM_each_sample(HaplotypePair prev_hap,
													size_t sample_index,
													bool exec) {
	const vector<int>	int_gts = get_int_gts(sample_index);
	const vector<Haplotype>	haps_mat = collect_haplotypes_mat(sample_index);
	const vector<Haplotype>	haps_pat = collect_haplotypes_pat(sample_index);
	// select a combination with seed
	const int	seed = get_record(0)->get_pos();
	const HaplotypePair	hap = Haplotype::impute(int_gts, haps_mat, haps_pat,
																prev_hap, seed);
	if(exec)	// actually impute?
		this->set_haplotype(hap, sample_index);
	return hap;
}

void VCFImputable::set_haplotype(HaplotypePair hap, size_t sample_index) {
	const Haplotype&	hap_mat = hap.first;
	const Haplotype&	hap_pat = hap.second;
	for(size_t i = 0; i < hap_mat.hap.size(); ++i) {
		const int	gt = Genotype::from_alleles(hap_mat.hap[i], hap_pat.hap[i]);
		get_record(i)->set_geno(sample_index, gt);
	}
}
