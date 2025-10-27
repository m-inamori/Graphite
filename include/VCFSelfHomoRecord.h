#ifndef __VCFHOMOHOMO
#define __VCFHOMOHOMO

#include "VCFImpSelfRecord.h"


//////////////////// VCFSelfHomoRecord ////////////////////

class VCFSelfHomoRecord : public VCFImpSelfRecord {
public:
	VCFSelfHomoRecord(ll pos, const std::vector<int>& geno,
						int i, WrongType type, ParentComb c) :
								VCFImpSelfRecord(pos, geno, i, type, c) { }
	~VCFSelfHomoRecord() { }
	
	bool is_imputable() const override { return false; }
	SelfFillType get_fill_type() const override { return SelfFillType::FILLED; }
	
	std::vector<int> gts() const;
	
	VCFSelfHomoRecord *impute() const;
};
#endif
