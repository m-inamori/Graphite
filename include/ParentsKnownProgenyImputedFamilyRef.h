#ifndef __PARENTSKNOWNPROGENYIMPUTEDFAMILYREF
#define __PARENTSKNOWNPROGENYIMPUTEDFAMILYREF

#include "VCF.h"
#include "VCFImpFamilyRecord.h"

class VCFImputable;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ParentsKnownProgenyImputedFamilyRef ////////////////////

namespace ParentsKnownProgenyImputedFamilyRef {
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
	VCFGeno *impute(const VCFSmall *orig_vcf,
							const VCFGeno *phased_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const std::vector<STRVEC>& imputed_progenies,
							const OptionSmall& op);
};
#endif
