#ifndef __ORPHAN
#define __ORPHAN

#include <vector>
#include "VCF.h"

class GenoRecord;
class VCFGenoBase;
class VCFOrphan;
class VCFOrphanRough;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// Orphan ////////////////////

namespace Orphan {
	VCFGenoBase *impute_samples(const STRVEC& samples,
							const std::vector<GenoRecord *>& records,
							const std::vector<std::vector<int>>& ref_haps,
							const VCFSmall *vcf, const OptionSmall& op);
	VCFGenoBase *impute(const STRVEC& samples,
						const VCFSmall *orig_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						const OptionSmall& op);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
													const OptionSmall& op);
	// upper NH which passes is_small
	std::size_t compute_upper_NH(std::size_t M, const OptionSmall& op);
};
#endif
