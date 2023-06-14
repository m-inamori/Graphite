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
	
	bool is_homohomo() const { return true; }
	bool is_imputable() const;
	FillType get_fill_type() const { 
		return is_imputable() ? FillType::IMPUTABLE : FillType::FILLED;
	}
	
	std::vector<std::string> gts() const;
	
	VCFHomoHomoRecord *impute() const;
};


//////////////////// VCFHomoHomo ////////////////////

class VCFHomoHomo: public VCFFamilyBase {
	const std::vector<VCFHomoHomoRecord *>	records;
	
public:
	VCFHomoHomo(const std::vector<STRVEC>& h, const STRVEC& s,
								std::vector<VCFHomoHomoRecord *> rs);
	
	
	///// virtual methods /////
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	std::vector<VCFHomoHomo *> impute() const;
};
#endif
