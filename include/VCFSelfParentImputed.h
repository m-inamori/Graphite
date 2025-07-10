#ifndef __VCFSELFPARENTIMPUTED
#define __VCFSELFPARENTIMPUTED

#include <vector>
#include "VCF.h"

class Map;
class SelfProgenyImputer;


//////////////////// VCFSelfParentImputed ////////////////////

class VCFSelfParentImputed : public VCFBase, public VCFSmallBase {
	const std::vector<VCFRecord *>	records;
	std::vector<SelfProgenyImputer *>	imputers;
	
public:
	VCFSelfParentImputed(const std::vector<STRVEC>& header, const STRVEC& s,
									const std::vector<VCFRecord *>& records,
									const Map& map_, double w);
	VCFSelfParentImputed(const VCFSelfParentImputed&) = delete;
	VCFSelfParentImputed& operator=(const VCFSelfParentImputed&) = delete;
	~VCFSelfParentImputed();
	
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
