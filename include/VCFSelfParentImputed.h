#ifndef __VCFSELFPARENTIMPUTED
#define __VCFSELFPARENTIMPUTED

#include <vector>
#include "SelfProgenyImputer.h"

class Map;


//////////////////// VCFSelfParentImputed ////////////////////

class VCFSelfParentImputed : public VCFGenoBase {
	const std::vector<GenoRecord *>	records;
	SelfProgenyImputer	prog_imputer;
	
public:
	VCFSelfParentImputed(const STRVEC& s,
							const std::vector<GenoRecord *>& records,
							const Map& map_, double w, const VCFSmall *vcf);
	VCFSelfParentImputed(const VCFSelfParentImputed&) = delete;
	VCFSelfParentImputed& operator=(const VCFSelfParentImputed&) = delete;
	~VCFSelfParentImputed();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
	
private:
	std::size_t num_progenies() const { return get_samples().size() - 1; }
};
#endif
