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

class VCFHeteroHeteroLite: public VCFFamilyBase {
	const std::vector<VCFHeteroHeteroLiteRecord *>	records;
	
public:
	VCFHeteroHeteroLite(const std::vector<STRVEC>& h, const STRVEC& s,
								std::vector<VCFHeteroHeteroLiteRecord *> rs);
	~VCFHeteroHeteroLite() { }
	
	///// virtual methods /////
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
};
#endif
