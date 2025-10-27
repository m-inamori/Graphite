#ifndef __VCFHETEROHETEROLITE
#define __VCFHETEROHETEROLITE

#include "VCFImpFamilyRecord.h"


//////////////////// VCFHeteroHeteroLiteRecord ////////////////////

class VCFHeteroHeteroLiteRecord : public VCFImpFamilyRecord {
public:
	VCFHeteroHeteroLiteRecord(ll pos, const std::vector<int>& geno,
									int i, WrongType type, ParentComb c) :
								VCFImpFamilyRecord(pos, geno, i, type, c) { }
	~VCFHeteroHeteroLiteRecord() { }
	
	bool is_homohomo() const override { return false; }
	bool is_imputable() const override { return false; }
	FillType get_fill_type() const override { return FillType::IMPUTABLE; }
};


//////////////////// VCFHeteroHeteroLite ////////////////////

class VCFHeteroHeteroLite: public VCFFamilyBase {
	const std::vector<VCFHeteroHeteroLiteRecord *>	records;
	
public:
	VCFHeteroHeteroLite(const STRVEC& s,
						const std::vector<VCFHeteroHeteroLiteRecord *>& rs,
						const VCFSmall *vcf);
	~VCFHeteroHeteroLite() { }
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
};
#endif
