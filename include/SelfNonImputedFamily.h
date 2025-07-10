#ifndef __SELFNONIMPUTEDFAMILY
#define __SELFNONIMPUTEDFAMILY

#include <vector>

class VCFSmallBase;
class VCFSmall;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// SelfNonImputedFamily ////////////////////

namespace SelfNonImputedFamily {
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
							const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	VCFSmallBase *impute(const VCFSmall *orig_vcf,
						 const VCFSmall *imputed_vcf,
						 const std::vector<std::vector<int>>& ref_haps,
						 const std::vector<const KnownFamily *>& families,
						 const OptionSmall& op);
};
#endif
