#ifndef __ONEKNOWNFAMILY
#define __ONEKNOWNFAMILY

#include <vector>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFImputable;
class VCFFamilyRecord;
class Family;
class KnownFamily;
class OptionSmall;


//////////////////// OneKnownFamily ////////////////////

namespace OneKnownFamily {
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
										std::size_t L, const OptionSmall& op);
	std::size_t compute_upper_NH(const Family *family, std::size_t M,
										std::size_t L, const OptionSmall& op);
	VCFImputable *create_family_vcf(const KnownFamily *family,
								const std::vector<VCFFamilyRecord *>& records,
								std::size_t num_families,
								const std::vector<std::vector<int>>& ref_haps,
								const VCFSmall *orig_vcf,
								const OptionSmall& op);
	VCFGenoBase *impute(const VCFSmall *orig_vcf, const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const OptionSmall& op);
};
#endif
