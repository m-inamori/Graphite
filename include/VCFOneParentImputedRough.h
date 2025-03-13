#ifndef __VCFONEPARENTIMPUTEDROUGH
#define __VCFONEPARENTIMPUTEDROUGH

#include "VCFFamily.h"
#include "VCFOneParentImputedBase.h"

class Map;
class ParentImputer;
class ProgenyImputer;


//////////////////// VCFOneParentImputedRough ////////////////////

class VCFOneParentImputedRough : public VCFOneParentImputedBase {
	const std::vector<VCFFamilyRecord *>	ref_records;
	ParentImputer	*parent_imputer;
	ProgenyImputer	*prog_imputer;
	
public:
	VCFOneParentImputedRough(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_imputed, const Map& map_, double w);
	VCFOneParentImputedRough(const VCFOneParentImputedRough&) = delete;
	VCFOneParentImputedRough& operator=(const VCFOneParentImputedRough&) = delete;
	~VCFOneParentImputedRough();
	
	///// virtual methods /////
	void impute() override;
	std::size_t amount() const override;
};
#endif
