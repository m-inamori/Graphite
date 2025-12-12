#ifndef __VCFONEPARENTIMPUTEDROUGH
#define __VCFONEPARENTIMPUTEDROUGH

#include "VCFFamily.h"
#include "VCFImputable.h"
#include "ParentImputer.h"
#include "ProgenyImputer.h"

class Map;


//////////////////// VCFOneParentImputedRough ////////////////////

class VCFOneParentImputedRough : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	ParentImputer	parent_imputer;
	ProgenyImputer	prog_imputer;
	
public:
	VCFOneParentImputedRough(const STRVEC& s,
								const std::vector<VCFFamilyRecord *>& records,
								const std::vector<std::vector<int>>& ref_haps,
								bool is_mat_imputed, const Map& map_, double w,
								const VCFSmall *vcf);
	VCFOneParentImputedRough(const VCFOneParentImputedRough&) = delete;
	VCFOneParentImputedRough& operator=(const VCFOneParentImputedRough&) = delete;
	~VCFOneParentImputedRough();
	
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
	std::size_t amount() const override;
	void impute() override;
};
#endif
