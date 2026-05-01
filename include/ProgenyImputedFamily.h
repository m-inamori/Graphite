#ifndef __PROJECTIMPUTEDFAMILY
#define __PROJECTIMPUTEDFAMILY

#include <vector>
#include <string>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFProgenyImputed;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ProgenyImputedFamily ////////////////////

namespace ProgenyImputedFamily {
	std::vector<std::vector<std::string>> collect_samples(
				const std::vector<const KnownFamily *>& families,
				const std::vector<std::vector<std::string>>& imputed_progenies);
	VCFGeno *impute(const VCFSmall *orig_vcf, const VCFGenoBase *imputed_vcf,
				const std::vector<const KnownFamily *>& families,
				const std::vector<std::vector<std::string>>& imputed_progenies,
				const std::vector<std::vector<int>>& ref_haps,
				const OptionSmall& op);
};
#endif
