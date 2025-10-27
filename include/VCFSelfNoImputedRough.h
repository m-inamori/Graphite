#ifndef __VCFSELFNOIMPUTEDROUGH
#define __VCFSELFNOIMPUTEDROUGH

#include <vector>
#include "VCFGeno.h"
#include "SelfParentImputerLessImputed.h"

class VCFSmall;
class GenoRecord;
class Map;
class SelfProgenyImputer;


//////////////////// VCFSelfNoImputedRough ////////////////////

class VCFSelfNoImputedRough : public VCFGenoBase, VCFMeasurable {
	const std::vector<GenoRecord *>	records;
	SelfParentImputerLessImputed	imputer;
	std::vector<SelfProgenyImputer *>	imputers;
	
public:
	VCFSelfNoImputedRough(const STRVEC& s,
							const std::vector<GenoRecord *>& records,
							const std::vector<std::vector<int>>& ref_hs,
							const Map& map_, double w, const VCFSmall *vcf);
	VCFSelfNoImputedRough(const VCFSelfNoImputedRough&) = delete;
	VCFSelfNoImputedRough& operator=(const VCFSelfNoImputedRough&) = delete;
	~VCFSelfNoImputedRough();
	
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
