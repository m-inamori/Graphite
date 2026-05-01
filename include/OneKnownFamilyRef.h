#ifndef __ONEKNOWNFAMILYREF
#define __ONEKNOWNFAMILYREF

#include <vector>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFOneParentKnown;
class Family;
class KnownFamily;
class OptionSmall;


//////////////////// OneKnownFamilyRef ////////////////////

namespace OneKnownFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf, const VCFGeno *phased_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						const std::vector<const KnownFamily *>& families,
						const OptionSmall& op);
};
#endif
