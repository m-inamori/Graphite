#ifndef __VCFSELFNOIMPUTEDROUGH
#define __VCFSELFNOIMPUTEDROUGH

#include <vector>
#include "VCFGeno.h"
#include "VCFSelfImputable.h"
#include "SelfParentImputerLessImputed.h"
#include "SelfProgenyImputer.h"

class VCFSmall;
class GenoRecord;
class Map;


//////////////////// VCFSelfNoImputedRough ////////////////////

class VCFSelfNoImputedRough : public VCFSelfImputable {
	const std::vector<GenoRecord *>	records;
	SelfParentImputerLessImputed	parent_imputer;
	SelfProgenyImputer				prog_imputer;
	
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
	
	///// virtual methods for VCFSelfImputable /////
	std::size_t amount() const override { return 1; }
	void impute() override;
	
private:
	std::size_t num_progenies() const { return get_samples().size() - 1; }
};
#endif
