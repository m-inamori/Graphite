#ifndef __ONEKNOWNFAMILY
#define __ONEKNOWNFAMILY

#include <vector>
#include "VCF.h"

class VCFOneParentKnown;
class Family;
class KnownFamily;
class Map;


//////////////////// OneKnownFamily ////////////////////

namespace OneKnownFamily {
	struct ConfigThread {
		int	first;
		int	num_threads;
		const std::vector<VCFOneParentKnown *>&	vcfs;
		
		ConfigThread(int i, int n,
						const std::vector<VCFOneParentKnown *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	VCFSmallBase *impute(const VCFSmall *orig_vcf, const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps, std::size_t L);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFOneParentKnown *>& vcfs, int T);
	VCFSmallBase *impute(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
};
#endif
