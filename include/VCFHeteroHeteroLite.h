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
	
	bool is_homohomo() const override { return false; }
	bool is_imputable() const override { return false; }
	FillType get_fill_type() const override { return FillType::IMPUTABLE; }
};


//////////////////// VCFHeteroHeteroLite ////////////////////

class VCFHeteroHeteroLite: public VCFBase, public VCFSmallBase,
											public VCFFamilyBase {
	const std::vector<VCFHeteroHeteroLiteRecord *>	records;
	
public:
	VCFHeteroHeteroLite(const std::vector<STRVEC>& h, const STRVEC& s,
						const std::vector<VCFHeteroHeteroLiteRecord *>& rs);
	~VCFHeteroHeteroLite() { }
	
	///// virtual methods /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
};
#endif
