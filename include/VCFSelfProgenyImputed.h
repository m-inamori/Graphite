#ifndef __VCFSELFPROGENYIMPUTED
#define __VCFSELFPROGENYIMPUTED

#include <vector>
#include "VCF.h"

class Map;
class SelfParentImputer;
class SelfProgenyImputer;


//////////////////// VCFSelfProgenyImputed ////////////////////

class VCFSelfProgenyImputed : public VCFBase, public VCFSmallBase {
	const std::vector<VCFRecord *>	records;
	std::size_t	ic;
	SelfParentImputer	*imputer;
	std::vector<SelfProgenyImputer *>	imputers;
	
public:
	VCFSelfProgenyImputed(const std::vector<STRVEC>& header, const STRVEC& s,
								const std::vector<VCFRecord *>& records,
								const std::vector<std::vector<int>>& ref_hs,
								std::size_t ic_, const Map& map_, double w);
	VCFSelfProgenyImputed(const VCFSelfProgenyImputed&) = delete;
	VCFSelfProgenyImputed& operator=(const VCFSelfProgenyImputed&) = delete;
	~VCFSelfProgenyImputed();
	
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
