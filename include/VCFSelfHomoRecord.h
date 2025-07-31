#ifndef __VCFHOMOHOMO
#define __VCFHOMOHOMO

#include "VCFImpSelfRecord.h"


//////////////////// VCFSelfHomoRecord ////////////////////

class VCFSelfHomoRecord : public VCFImpSelfRecord {
public:
	VCFSelfHomoRecord(const STRVEC& v, const STRVEC& samples,
						int i, WrongType type, ParentComb c) :
								VCFImpSelfRecord(v, samples, i, type, c) { }
	~VCFSelfHomoRecord() { }
	
	bool is_imputable() const override { return false; }
	SelfFillType get_fill_type() const override { return SelfFillType::FILLED; }
	
	std::vector<std::string> gts() const;
	
	VCFSelfHomoRecord *impute() const;
};
#endif
