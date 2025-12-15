#ifndef __VCFORPHANROUGH
#define __VCFORPHANROUGH

#include "VCFGeno.h"
#include "OrphanImputer.h"

class VCFSmall;
class Map;


//////////////////// VCFOrphanRough ////////////////////

class VCFOrphanRough : public VCFGenoBase {
public:
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		VCFOrphanRough	*vcf;
		
		ConfigThread(std::size_t i, std::size_t n, VCFOrphanRough *vcf_) :
									first(i), num_threads(n), vcf(vcf_) { }
	};
	
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
	void impute_each(std::size_t i);
	void impute(int num_threads);
	
private:
	std::vector<OrphanImputer *> create_imputers(
			const std::vector<std::vector<std::vector<int>>>& ref_haps_table,
			const Map& map_, double w);
	
public:
	static void impute_small_in_thread(void *config);
};
#endif
