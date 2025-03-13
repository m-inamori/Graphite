#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFOneParentImputedBase.h"
#include "VCFFamily.h"
#include "Map.h"

class ParentImputer;
class ProgenyImputer;


//////////////////// VCFOneParentKnown ////////////////////

class VCFOneParentKnown : public VCFBase, public VCFSmallBase,
											public VCFFamilyBase {
	const std::vector<VCFFamilyRecord *>	records;
	const bool is_mat_known;
	ParentImputer	*parent1_imputer;
	ParentImputer	*parent2_imputer;
	ProgenyImputer	*prog_imputer;
	
public:
	VCFOneParentKnown(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_known, const Map& map_, double w);
	VCFOneParentKnown(const VCFOneParentKnown&) = delete;
	VCFOneParentKnown& operator=(const VCFOneParentKnown&) = delete;
	~VCFOneParentKnown();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
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
	const std::string& get_known_parent() const;
	
	void impute_known_parent();
	void impute();
};
#endif
