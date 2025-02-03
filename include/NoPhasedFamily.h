#ifndef __NOPHASEDFAMILY
#define __NOPHASEDFAMILY

#include <vector>
#include "VCF.h"

class VCFFamily;
class VCFNoParentImputed;
class Family;
class KnownFamily;
class Map;


//////////////////// NoPhasedFamily ////////////////////

namespace NoPhasedFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFNoParentImputed *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFNoParentImputed *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	void create_in_thread(void *config);
	
	VCFSmallBase *impute(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFNoParentImputed *>& vcfs, int T);
};
#endif
