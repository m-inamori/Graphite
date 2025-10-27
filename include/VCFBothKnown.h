#ifndef __VCFNOPARENTIMPUTED
#define __VCFNOPARENTIMPUTED

#include "VCFFamily.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"

class Map;
class ParentImputer;
class ProgenyImputer;


//////////////////// VCFBothKnown ////////////////////

class VCFBothKnown : public VCFFamilyBase {
private:
	std::vector<VCFFamilyRecord *>	records;
	ParentImputer	mat_imputer;
	ParentImputer	pat_imputer;
	ProgenyImputer	prog_imputer;
	
public:
	VCFBothKnown(const STRVEC& s,
					const std::vector<VCFFamilyRecord *>& records,
					const std::vector<std::vector<int>>& ref_haps,
					const Map& map_, double w, const VCFSmall *vcf);
	VCFBothKnown(const VCFBothKnown&) = delete;
	VCFBothKnown& operator=(const VCFBothKnown&) = delete;
	~VCFBothKnown();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
	void clear_records() { records.clear(); }
};
#endif
