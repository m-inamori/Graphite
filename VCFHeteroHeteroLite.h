#ifndef __VCFHETEROHETEROLITE
#define __VCFHETEROHETEROLITE

#include "VCFImpFamily.h"


//////////////////// VCFHeteroHeteroLiteRecord ////////////////////

class VCFHeteroHeteroLiteRecord : public VCFImpFamilyRecord {
public:
	VCFHeteroHeteroLiteRecord(const STRVEC& v, const STRVEC& samples,
									int i, WrongType type, ParentComb c) :
								VCFImpFamilyRecord(v, samples, i, type, c) { }
	~VCFHeteroHeteroLiteRecord() { }
	
	bool is_homohomo() const { return false; }
	bool is_imputable() const { return false; }
	FillType get_fill_type() const { return FillType::IMPUTABLE; }
};


//////////////////// VCFHeteroHeteroLite ////////////////////

class VCFHeteroHeteroLite: public VCFFamily {
	const std::vector<VCFHeteroHeteroLiteRecord *>	hehe_records;
	
public:
	VCFHeteroHeteroLite(const std::vector<STRVEC>& h, const STRVEC& s,
								std::vector<VCFHeteroHeteroLiteRecord *> rs);
};
#endif
