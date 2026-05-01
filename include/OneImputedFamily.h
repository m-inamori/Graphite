#ifndef __ONEIMPUTEDFAMILY
#define __ONEIMPUTEDFAMILY

#include <vector>

class VCFSmallBase;
class VCFSmall;
class VCFGeno;
class VCFImputable;
class VCFFamilyRecord;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// OneImputedFamily ////////////////////

namespace OneImputedFamily {
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
						const std::vector<std::vector<int>>& ref_haps,
						std::size_t L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
									std::size_t L, const OptionSmall& op);
	VCFImputable *create_family_vcf(const KnownFamily *family,
								bool is_mat_imputed,
								const std::vector<VCFFamilyRecord *>& records,
								std::size_t num_families,
								const std::vector<std::vector<int>>& ref_haps,
								const VCFSmall *vcf,
								const OptionSmall& op);
	VCFGenoBase *impute(const VCFSmall *orig_vcf, const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const OptionSmall& op);
};
#endif
