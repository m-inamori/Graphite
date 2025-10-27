#ifndef __VCFNOPARENTKNOWN
#define __VCFNOPARENTKNOWN

#include "VCFFamily.h"
#include "ParentImputer.h"
#include "ProgenyImputer.h"

class Map;


//////////////////// VCFNoParentKnown ////////////////////

class VCFNoParentKnown : public VCFFamilyBase {
private:
	std::vector<VCFFamilyRecord *>	records;
	ParentImputer	mat_imputer;
	ParentImputer	pat_imputer;
	ProgenyImputer	prog_imputer;
	
public:
	VCFNoParentKnown(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w, const VCFSmall *vcf);
	VCFNoParentKnown(const VCFNoParentKnown&) = delete;
	VCFNoParentKnown& operator=(const VCFNoParentKnown&) = delete;
	~VCFNoParentKnown();
	
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
