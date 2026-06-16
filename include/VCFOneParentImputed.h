#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFImputable.h"
#include "ImputerByParentProgeny.h"
#include "Map.h"


//////////////////// VCFOneParentImputed ////////////////////

class VCFOneParentImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	ImputerByParentProgeny	imputer;
	const bool	should_impute_mat;
	
public:
	VCFOneParentImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool should_impute_mat_, const Map& map_,
						double w, const VCFSmall *vcf);
	VCFOneParentImputed(const VCFOneParentImputed&) = delete;
	VCFOneParentImputed& operator=(const VCFOneParentImputed&) = delete;
	~VCFOneParentImputed();
	
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
