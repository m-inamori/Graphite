#ifndef __VCFORPHANROUGH
#define __VCFORPHANROUGH

#include "VCFGeno.h"
#include "OrphanImputer.h"

class VCFSmall;
class Map;


//////////////////// VCFOrphanRough ////////////////////

class VCFOrphanRough : public VCFGenoBase {
private:
	std::vector<GenoRecord *>	records;
	std::vector<OrphanImputer *>	imputers;
	
public:
	VCFOrphanRough(const STRVEC& s, const std::vector<GenoRecord *>& records,
			const std::vector<std::vector<std::vector<int>>>& ref_haps_table,
			const Map& map_, double w, const VCFSmall *vcf);
	VCFOrphanRough(const VCFOrphanRough&) = delete;
	VCFOrphanRough& operator=(const VCFOrphanRough&) = delete;
	~VCFOrphanRough();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute(std::size_t io);
	
private:
	std::vector<OrphanImputer *> create_imputers(
			const std::vector<std::vector<std::vector<int>>>& ref_haps_table,
			const Map& map_, double w);
};
#endif
