#ifndef __VCFHOMOHOMO
#define __VCFHOMOHOMO

#include "VCFImpFamilyRecord.h"


//////////////////// VCFHomoHomoRecord ////////////////////

class VCFHomoHomoRecord : public VCFImpFamilyRecord {
public:
	VCFHomoHomoRecord(ll pos, const std::vector<int>& geno,
						int i, WrongType type, ParentComb c) :
								VCFImpFamilyRecord(pos, geno, i, type, c) { }
	~VCFHomoHomoRecord() { }
	
	bool is_homohomo() const override { return true; }
	bool is_imputable() const override {
		return comb == ParentComb::P00x11 && is_mat_NA();
	}
	FillType get_fill_type() const override { 
		return is_imputable() ? FillType::IMPUTABLE : FillType::FILLED;
	}
	
	std::vector<int> gts() const;
	
	VCFHomoHomoRecord *impute() const;
};


//////////////////// VCFHomoHomo ////////////////////

class VCFHomoHomo: public VCFFamilyBase {
	const std::vector<VCFHomoHomoRecord *>	records;
	
public:
	VCFHomoHomo(const STRVEC& s, const std::vector<VCFHomoHomoRecord *>& rs,
														const VCFSmall *vcf);
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override { return records[i]; }
	
	///// virtual methods for VCFImputable /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	std::vector<VCFHomoHomo *> impute() const;
};
#endif
