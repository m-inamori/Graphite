#ifndef __VCFJUNKRECORD
#define __VCFJUNKRECORD

#include "VCFImpFamilyRecord.h"


//////////////////// VCFJunkRecord ////////////////////

class VCFJunkRecord : public VCFImpFamilyRecord {
public:
	VCFJunkRecord(ll pos, const std::vector<int>& geno,
										int i, WrongType type) :
					VCFImpFamilyRecord(pos, geno, i, type, ParentComb::PNA) { }
	~VCFJunkRecord() { }
	
	bool is_homohomo() const override { return false; }
	bool is_imputable() const override { return false; }
	FillType get_fill_type() const override { return FillType::UNABLE; }
};
#endif
