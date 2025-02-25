#ifndef __VCFNOPARENTKNOWN
#define __VCFNOPARENTKNOWN

#include "VCFFamily.h"

class Map;
class ParentImputer;
class ProgenyImputer;


//////////////////// VCFNoParentKnown ////////////////////

class VCFNoParentKnown : public VCFBase, public VCFSmallBase {
private:
	std::vector<VCFFamilyRecord *>	records;
	ParentImputer	*mat_imputer;
	ParentImputer	*pat_imputer;
	ProgenyImputer	*prog_imputer;
	
public:
	VCFNoParentKnown(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w);
	~VCFNoParentKnown();
	
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
	void clear_records() { records.clear(); }
};
#endif
