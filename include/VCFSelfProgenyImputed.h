#ifndef __VCFSELFPROGENYIMPUTED
#define __VCFSELFPROGENYIMPUTED

#include <vector>
#include "VCF.h"
#include "SelfParentImputer.h"

class Map;
class SelfProgenyImputer;


//////////////////// VCFSelfProgenyImputed ////////////////////

class VCFSelfProgenyImputed : public VCFGenoBase {
	const std::vector<GenoRecord *>	records;
	std::size_t	ic;
	SelfParentImputer	imputer;
	std::vector<SelfProgenyImputer *>	imputers;
	
public:
	VCFSelfProgenyImputed(const STRVEC& s,
							const std::vector<GenoRecord *>& records,
							const std::vector<std::vector<int>>& ref_hs,
							std::size_t ic_, const Map& map_, double w,
							const VCFSmall *vcf);
	VCFSelfProgenyImputed(const VCFSelfProgenyImputed&) = delete;
	VCFSelfProgenyImputed& operator=(const VCFSelfProgenyImputed&) = delete;
	~VCFSelfProgenyImputed();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
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
