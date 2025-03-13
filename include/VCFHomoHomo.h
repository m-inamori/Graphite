#ifndef __VCFHOMOHOMO
#define __VCFHOMOHOMO

#include "VCFImpFamily.h"


//////////////////// VCFHomoHomoRecord ////////////////////

class VCFHomoHomoRecord : public VCFImpFamilyRecord {
public:
	VCFHomoHomoRecord(const STRVEC& v, const STRVEC& samples,
						int i, WrongType type, ParentComb c) :
								VCFImpFamilyRecord(v, samples, i, type, c) { }
	~VCFHomoHomoRecord() { }
	
	bool is_homohomo() const override { return true; }
	bool is_imputable() const override;
	FillType get_fill_type() const override { 
		return is_imputable() ? FillType::IMPUTABLE : FillType::FILLED;
	}
	
	std::vector<std::string> gts() const;
	
	VCFHomoHomoRecord *impute() const;
};


//////////////////// VCFHomoHomo ////////////////////

class VCFHomoHomo: public VCFBase, public VCFFamilyBase {
	const std::vector<VCFHomoHomoRecord *>	records;
	
public:
	VCFHomoHomo(const std::vector<STRVEC>& h, const STRVEC& s,
							const std::vector<VCFHomoHomoRecord *>& rs);
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	
	///// virtual methods for VCFImputable /////
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	std::vector<VCFHomoHomo *> impute() const;
};
#endif
