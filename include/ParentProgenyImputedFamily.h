#ifndef __PARENTPROGENYIMPUTEDFAMILY
#define __PARENTPROGENYIMPUTEDFAMILY

#include <vector>
#include <string>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFProgenyImputed;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ParentProgenyImputedFamily ////////////////////

namespace ParentProgenyImputedFamily {
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGenoBase *imputed_vcf,
					const std::vector<const KnownFamily *>& families,
					const std::vector<std::string>& non_imputed_parents,
					const std::vector<std::vector<int>>& ref_haps,
					const OptionSmall& op);
};
#endif
