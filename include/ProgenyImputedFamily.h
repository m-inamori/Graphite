#ifndef __PROJECTIMPUTEDFAMILY
#define __PROJECTIMPUTEDFAMILY

#include <vector>
#include <string>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFProgenyImputed;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ProgenyImputedFamily ////////////////////

namespace ProgenyImputedFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFProgenyImputed *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFProgenyImputed *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFProgenyImputed *>& v, int T);
	VCFGeno *impute(const VCFSmall *orig_vcf, const VCFGenoBase *imputed_vcf,
				const std::vector<const KnownFamily *>& families,
				const std::vector<std::vector<std::string>>& imputed_progenies,
				const std::vector<std::vector<int>>& ref_haps,
				const OptionSmall& op);
};
#endif
