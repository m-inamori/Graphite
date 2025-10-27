#include "../include/common.h"
#include "../include/VCFHomoHomo.h"

using namespace std;


//////////////////// VCFHomoHomoRecord ////////////////////

vector<int> VCFHomoHomoRecord::gts() const {
	if(this->comb == ParentComb::P00x00)
		return vector<int>(geno.size(), Genotype::PH_00);
	else if(this->comb == ParentComb::P11x11)
		return vector<int>(geno.size(), Genotype::PH_11);
	
	// 0/0 x 1/1 is difficult
	if((this->is_mat_ref_homo() && !this->is_pat_ref_homo()) ||
					(!this->is_mat_alt_homo() && this->is_pat_alt_homo())) {
		vector<int>	gts(geno.size(), Genotype::PH_01);
		gts[0] = Genotype::PH_00;
		gts[1] = Genotype::PH_11;
		return gts;
	}
	else if((this->is_mat_alt_homo() && !this->is_pat_alt_homo()) ||
					(!this->is_mat_ref_homo() && this->is_pat_ref_homo())) {
		vector<int>	gts(geno.size(), Genotype::PH_10);
		gts[0] = Genotype::PH_11;
		gts[1] = Genotype::PH_00;
		return gts;
	}
	else {
		vector<int>	gts(geno.size(), Genotype::UN_01);
		gts[0] = Genotype::NA;
		gts[1] = Genotype::NA;
		return gts;
	}
}

VCFHomoHomoRecord *VCFHomoHomoRecord::impute() const {
	return new VCFHomoHomoRecord(this->pos, this->gts(), this->index,
												this->wrong_type, this->comb);
}


//////////////////// VCFHomoHomo ////////////////////

VCFHomoHomo::VCFHomoHomo(const STRVEC& s, const vector<VCFHomoHomoRecord *>& rs,
														const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf), records(rs) { }

vector<VCFHomoHomo *> VCFHomoHomo::impute() const {
	vector<VCFHomoHomoRecord *>	new_records;
	for(auto p = records.begin(); p != records.end(); ++p)
		new_records.push_back((*p)->impute());
	// HeteroHetero returns a vector and matches it
	auto	*vcf = new VCFHomoHomo(this->get_samples(),
									new_records, this->get_ref_vcf());
	return vector<VCFHomoHomo *>(1U, vcf);
}
