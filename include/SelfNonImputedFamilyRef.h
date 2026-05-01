#ifndef __SELFNONIMPUTEDFAMILYREF
#define __SELFNONIMPUTEDFAMILYREF

#include <vector>

class VCFSmall;
class VCFGeno;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// SelfNonImputedFamilyRef ////////////////////

namespace SelfNonImputedFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGeno *imputed_vcf,
					const std::vector<std::vector<int>>& ref_haps,
					const std::vector<const KnownFamily *>& families,
					const OptionSmall& op);
};
#endif
