#ifndef __VCFORPHAN
#define __VCFORPHAN

#include "VCFGeno.h"
#include "OrphanImputer.h"

class VCFSmall;
class Map;


//////////////////// VCFOrphan ////////////////////

class VCFOrphan : public VCFGenoBase {
private:
	std::vector<GenoRecord *>	records;
	OrphanImputer	imputer;
	
public:
	VCFOrphan(const STRVEC& s, const std::vector<GenoRecord *>& records,
				const std::vector<std::vector<int>>& ref_haps,
				const Map& map_, double w, const VCFSmall *vcf);
	VCFOrphan(const VCFOrphan&) = delete;
	VCFOrphan& operator=(const VCFOrphan&) = delete;
	~VCFOrphan();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute(std::size_t io);
};
#endif
