#ifndef __VCFJUNKRECORD
#define __VCFJUNKRECORD

#include "VCFImpFamily.h"


//////////////////// VCFJunkRecord ////////////////////

class VCFJunkRecord : public VCFImpFamilyRecord {
public:
	VCFJunkRecord(const STRVEC& v, const STRVEC& samples,
										int i, WrongType type) :
					VCFImpFamilyRecord(v, samples, i, type, ParentComb::PNA) { }
	~VCFJunkRecord() { }
	
	bool is_homohomo() const { return false; }
	bool is_imputable() const { return false; }
	FillType get_fill_type() const { return FillType::UNABLE; }
};
#endif
