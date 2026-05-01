#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFImputable.h"
#include "ParentImputer.h"
#include "ProgenyImputer.h"
#include "Map.h"


//////////////////// VCFOneParentKnown ////////////////////

class VCFOneParentKnown :  public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const bool is_mat_known;
	ParentImputer	parent_imputer;
	
public:
	VCFOneParentKnown(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_known, const Map& map_, double w,
						const VCFSmall *vcf);
	VCFOneParentKnown(const VCFOneParentKnown&) = delete;
	VCFOneParentKnown& operator=(const VCFOneParentKnown&) = delete;
	~VCFOneParentKnown();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	void impute() override;
	std::size_t amount() const override { return 1; }
};
#endif
