#ifndef __ONEIMPUTEDFAMILYREF
#define __ONEIMPUTEDFAMILYREF

#include <vector>

class VCFSmallBase;
class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFImputable;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// OneImputedFamilyRef ////////////////////

namespace OneImputedFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf, const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const OptionSmall& op);
};
#endif
