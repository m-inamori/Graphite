#ifndef __SELFNONIMPUTEDFAMILY
#define __SELFNONIMPUTEDFAMILY

#include <vector>

class VCFSmall;
class VCFGeno;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// SelfNonImputedFamily ////////////////////

namespace SelfNonImputedFamily {
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
							const std::vector<std::vector<int>>& ref_haps,
							std::size_t L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
							std::size_t L, const OptionSmall& op);
	// upper NH which passes is_small_ref
	std::size_t compute_upper_NH(const Family *family, std::size_t M,
										std::size_t L, const OptionSmall& op);
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGeno *imputed_vcf,
					const std::vector<std::vector<int>>& ref_haps,
					const std::vector<const KnownFamily *>& families,
					const OptionSmall& op);
};
#endif
