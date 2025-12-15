#ifndef __VCFORPHAN
#define __VCFORPHAN

#include "VCFGeno.h"
#include "OrphanImputer.h"

class VCFSmall;
class Map;


//////////////////// VCFOrphan ////////////////////

class VCFOrphan : public VCFGenoBase {
public:
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		VCFOrphan	*vcf;
		
		ConfigThread(std::size_t i, std::size_t n, VCFOrphan *vcf_) :
									first(i), num_threads(n), vcf(vcf_) { }
	};
	
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
	void impute_each(std::size_t i);
	void impute(int num_threads);
	
public:
	static void impute_small_in_thread(void *config);
};
#endif
