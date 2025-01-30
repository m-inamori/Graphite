#ifndef __VCFONEPARENTIMPUTEDROUGH
#define __VCFONEPARENTIMPUTEDROUGH

#include "VCFFamily.h"

class Map;
class ParentImputer;
class ProgenyImputer;


//////////////////// VCFOneParentImputedRough ////////////////////

class VCFOneParentImputedRough : public VCFBase, public VCFSmallBase {
	const std::vector<VCFFamilyRecord *>	records;
	ParentImputer	*parent_imputer;
	ProgenyImputer	*prog_imputer;
	
public:
	VCFOneParentImputedRough(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_imputed, const Map& map_, double w);
	~VCFOneParentImputedRough();
	
	///// virtual methods for VCFSmallBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
