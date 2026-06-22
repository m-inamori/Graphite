#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFImputable.h"
#include "ImputerByParentProgeny.h"
#include "Map.h"


//////////////////// VCFOneParentProgenyImputed ////////////////////

class VCFOneParentProgenyImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	ImputerByParentProgeny	imputer;
	const bool	should_impute_mat;
	
public:
	VCFOneParentProgenyImputed(const STRVEC& s,
								const std::vector<VCFFamilyRecord *>& records,
								const std::vector<std::vector<int>>& ref_haps,
								bool should_impute_mat_, const Map& map_,
								double w, const VCFSmall *vcf);
	VCFOneParentProgenyImputed(const VCFOneParentProgenyImputed&) = delete;
	VCFOneParentProgenyImputed&
				operator=(const VCFOneParentProgenyImputed&) = delete;
	~VCFOneParentProgenyImputed();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	std::size_t amount() const override;
	STRVEC imputed_samples() const override;
	void impute() override;
};
#endif
