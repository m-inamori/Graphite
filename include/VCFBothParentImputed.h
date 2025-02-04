#ifndef __VCFBOTHPARENTIMPUTED
#define __VCFBOTHPARENTIMPUTED

#include "VCFFamily.h"

class Map;
class ProgenyImputer;


//////////////////// VCFBothParentImputed ////////////////////

class VCFBothParentImputed : public VCFBase, public VCFSmallBase {
	const std::vector<VCFFamilyRecord *>	records;
	ProgenyImputer	*prog_imputer;
	
public:
	VCFBothParentImputed(const std::vector<STRVEC>& header, const STRVEC& s,
								const std::vector<VCFFamilyRecord *>& records,
								const Map& map_, double w);
	~VCFBothParentImputed();
	
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
