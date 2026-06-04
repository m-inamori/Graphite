#ifndef __PARENTSKNOWNPROGENYIMPUTEDFAMILY
#define __PARENTSKNOWNPROGENYIMPUTEDFAMILY

#include "VCF.h"
#include "VCFImpFamilyRecord.h"

class VCFImputable;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ParentsKnownProgenyImputedFamily ////////////////////

namespace ParentsKnownProgenyImputedFamily {
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
										int L, const OptionSmall& op);
	std::size_t compute_upper_NH(std::size_t M, std::size_t L,
												const OptionSmall& op);
	VCFImputable *create_family_vcf(
							const STRVEC& samples,
							const std::vector<VCFFamilyRecord *>& records,
							int num_families,
							const std::vector<std::vector<int>>& ref_haps,
							const VCFSmall *vcf,
							const OptionSmall& op);
	VCFGenoBase *impute(const VCFSmall *orig_vcf,
							const VCFGenoBase *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const std::vector<STRVEC>& imputed_progenies,
							const std::vector<std::vector<int>>& ref_haps,
							const OptionSmall& op);
};
#endif
