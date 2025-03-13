#ifndef __BOTHIMPUTEDFAMILY
#define __BOTHIMPUTEDFAMILY

#include <vector>
#include "VCF.h"

class VCFFamily;
class VCFBothParentImputed;
class Family;
class KnownFamily;
class Map;


//////////////////// BothImputedFamily ////////////////////

namespace BothImputedFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFBothParentImputed *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFBothParentImputed *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	void impute_small_in_thread(void *config);
	
	VCFSmallBase *impute(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& gmap, int num_threads);
	void impute_VCFs(std::vector<VCFBothParentImputed *>& vcfs, int T);
};
#endif
