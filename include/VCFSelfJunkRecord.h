#ifndef __VCFSELFJUNKRECORD
#define __VCFSELFJUNKRECORD

#include "VCFImpSelfRecord.h"


//////////////////// VCFSelfJunkRecord ////////////////////

class VCFSelfJunkRecord : public VCFImpSelfRecord {
public:
	VCFSelfJunkRecord(const STRVEC& v, const STRVEC& samples,
										int i, WrongType type) :
					VCFImpSelfRecord(v, samples, i, type, ParentComb::PNA) { }
	~VCFSelfJunkRecord() { }
	
	bool is_imputable() const override { return false; }
	SelfFillType get_fill_type() const override { return SelfFillType::UNABLE; }
};
#endif
