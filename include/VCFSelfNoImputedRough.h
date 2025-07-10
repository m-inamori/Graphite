#ifndef __VCFSELFNOIMPUTEDROUGH
#define __VCFSELFNOIMPUTEDROUGH

#include <vector>
#include "VCF.h"

class Map;
class SelfParentImputerLessImputed;
class SelfProgenyImputer;


//////////////////// VCFSelfNoImputedRough ////////////////////

class VCFSelfNoImputedRough : public VCFBase, public VCFSmallBase {
	const std::vector<VCFRecord *>	records;
	SelfParentImputerLessImputed	*imputer;
	std::vector<SelfProgenyImputer *>	imputers;
	
public:
	VCFSelfNoImputedRough(const std::vector<STRVEC>& header, const STRVEC& s,
								const std::vector<VCFRecord *>& records,
								const std::vector<std::vector<int>>& ref_hs,
								const Map& map_, double w);
	VCFSelfNoImputedRough(const VCFSelfNoImputedRough&) = delete;
	VCFSelfNoImputedRough& operator=(const VCFSelfNoImputedRough&) = delete;
	~VCFSelfNoImputedRough();
	
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
	
private:
	std::size_t num_progenies() const { return get_samples().size() - 1; }
	std::vector<SelfProgenyImputer *> create_imputers(
												const Map& map_, double w);
};
#endif
