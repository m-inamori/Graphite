#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFImputable.h"
#include "ParentProgenyImputer.h"
#include "Map.h"


//////////////////// VCFOneParentImputed ////////////////////

class VCFOneParentImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	ParentProgenyImputer	imputer;
	
public:
	VCFOneParentImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_imputed, const Map& map_,
						double w, const VCFSmall *vcf);
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
	void impute() override;
};
#endif
