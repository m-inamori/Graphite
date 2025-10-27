#ifndef __NOPHASEDFAMILY
#define __NOPHASEDFAMILY

#include <vector>
#include "VCF.h"

class VCFGenoBase;
class VCFGeno;
class VCFBothKnown;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// BothKnownFamily ////////////////////

namespace BothKnownFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFBothKnown *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFBothKnown *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	void create_in_thread(void *config);
	
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
											int L, const OptionSmall& op);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFBothKnown *>& vcfs, int T);
	VCFGenoBase *impute(const VCFSmall *orig_vcf,
						const VCFGeno *imputed_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						const std::vector<const KnownFamily *>& families,
						const OptionSmall& op);
};
#endif
