#ifndef __ONEIMPUTEDFAMILY
#define __ONEIMPUTEDFAMILY

#include <vector>
#include "VCF.h"

class VCFOneParentImputedBase;
class Family;
class KnownFamily;
class Map;


//////////////////// OneImputedFamily ////////////////////

namespace OneImputedFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFOneParentImputedBase *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFOneParentImputedBase *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	VCFSmallBase *impute(const VCFSmall *orig_vcf, const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
					std::size_t L);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
														std::size_t L);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFOneParentImputedBase *>& vcfs, int T);
	VCFSmallBase *impute(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
};
#endif
