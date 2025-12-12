#ifndef __VCFSELFNOIMPUTED
#define __VCFSELFNOIMPUTED

#include <vector>
#include "GenoRecord.h"
#include "SelfImputer.h"
#include "VCFGeno.h"
#include "Map.h"

class Map;
class SelfImputer;


//////////////////// VCFSelfNoImputed ////////////////////

class VCFSelfNoImputed : public VCFGenoBase, public VCFMeasurable {
	const std::vector<GenoRecord *>	records;
	SelfImputer	imputer;
	
public:
	VCFSelfNoImputed(const STRVEC& s, const std::vector<GenoRecord *>& records,
								const std::vector<std::vector<int>>& ref_hs,
								const Map& map_, double w, const VCFSmall *vcf);
	VCFSelfNoImputed(const VCFSelfNoImputed&) = delete;
	VCFSelfNoImputed& operator=(const VCFSelfNoImputed&) = delete;
	~VCFSelfNoImputed();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
