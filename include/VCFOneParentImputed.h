#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFOneParentImputedBase.h"
#include "ParentProgenyImputer.h"
#include "Map.h"


//////////////////// VCFOneParentImputed ////////////////////

class VCFOneParentImputed : public VCFOneParentImputedBase {
	const std::vector<std::vector<int>>&	ref_haps;
	ParentProgenyImputer	imputer;
	
public:
	VCFOneParentImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_imputed, const Map& map_,
						double w, const VCFSmall *vcf);
	~VCFOneParentImputed();
	
	///// virtual methods for VCFOneParentImputedBase /////
	std::size_t amount() const;
	void impute();
};
#endif
