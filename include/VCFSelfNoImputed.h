#ifndef __VCFSELFNOIMPUTED
#define __VCFSELFNOIMPUTED

#include <vector>
#include "VCF.h"

class Map;
class SelfImputer;


//////////////////// VCFSelfNoImputed ////////////////////

class VCFSelfNoImputed : public VCFBase, public VCFSmallBase {
	const std::vector<VCFRecord *>	records;
	SelfImputer	*imputer;
	
public:
	VCFSelfNoImputed(const std::vector<STRVEC>& header, const STRVEC& s,
								const std::vector<VCFRecord *>& records,
								const std::vector<std::vector<int>>& ref_hs,
								const Map& map_, double w);
	VCFSelfNoImputed(const VCFSelfNoImputed&) = delete;
	VCFSelfNoImputed& operator=(const VCFSelfNoImputed&) = delete;
	~VCFSelfNoImputed();
	
	///// virtual methods for VCFSmallBase /////
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
