#ifndef __VCFSELFHETERORECORD
#define __VCFSELFHETERORECORD

#include "VCFImpSelfRecord.h"


//////////////////// VCFSelfHeteroRecord ////////////////////

class VCFSelfHeteroRecord : public VCFImpSelfRecord {
public:
	VCFSelfHeteroRecord(ll pos, const std::vector<int>& geno,
							int i, WrongType type, ParentComb c) :
									VCFImpSelfRecord(pos, geno, i, type, c) { }
	
	std::size_t num_progenies() const { return num_samples() - 1; }
	
	bool is_imputable() const override {
		return this->wrong_type == WrongType::RIGHT;
	}
	SelfFillType get_fill_type() const override {
		return is_imputable() ? SelfFillType::P01 : SelfFillType::IMPUTABLE;
	}
	
	std::vector<int> progeny_gts() const;
	
	void set_haplo(int h);
	void set_int_gt_by_which_comes_from(int w1, int w2, std::size_t i);
};

#endif
