#ifndef __ORPHAN
#define __ORPHAN

#include <vector>
#include "VCF.h"

class VCFOrphan;
class Family;
class KnownFamily;
class Map;


//////////////////// Orphan ////////////////////

namespace Orphan {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		VCFOrphan	*vcf;
		
		ConfigThread(std::size_t i, std::size_t n, VCFOrphan *vcf_) :
									first(i), num_threads(n), vcf(vcf_) { }
	};
	
	void create_in_thread(void *config);
	void impute_small_in_thread(void *config);
	void impute_small_VCF(VCFOrphan *vcf, int T);
	
	VCFSmallBase *impute(const std::vector<std::string>& samples,
							const VCFSmall *orig_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const Map& gmap, int num_threads);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps);
};
#endif
