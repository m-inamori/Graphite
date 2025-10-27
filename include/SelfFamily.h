#ifndef __SELFFAMILY
#define __SELFFAMILY

#include <vector>
#include <string>

class VCFSmall;
class VCFGeno;
class KnownFamily;
class OptionSmall;


//////////////////// SelfFamily ////////////////////

namespace SelfFamily {
	VCFGeno *impute(const VCFSmall *orig_vcf,
						 const VCFGeno *imputed_vcf,
						 const std::vector<std::vector<int>>& ref_haps,
						 const std::vector<const KnownFamily *>& families,
						 const std::vector<std::string>& imputed_samples,
						 const OptionSmall& op);
};
#endif
