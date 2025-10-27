#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFFamily.h"
#include "ParentImputer.h"
#include "ProgenyImputer.h"
#include "Map.h"


//////////////////// VCFOneParentKnown ////////////////////

class VCFOneParentKnown :  public VCFFamilyBase {
	const std::vector<VCFFamilyRecord *>	records;
	const bool is_mat_known;
	ParentImputer	parent1_imputer;
	ParentImputer	parent2_imputer;
	ProgenyImputer	prog_imputer;
	
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
	
	///// non-virtual methods /////
	const std::string& get_known_parent() const;
	
	void impute_known_parent();
	void impute();
};
#endif
