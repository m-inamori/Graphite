#ifndef __ONEIMPUTEDFAMILY
#define __ONEIMPUTEDFAMILY

#include <vector>

class VCFSmallBase;
class VCFSmall;
class VCFGeno;
class VCFOneParentImputedBase;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


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
	
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
									std::size_t L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
									std::size_t L, const OptionSmall& op);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFOneParentImputedBase *>& vcfs, int T);
	VCFGenoBase *impute(const VCFSmall *orig_vcf, const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const OptionSmall& op);
};
#endif
