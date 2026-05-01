#ifndef __BOTHKNOWNFAMILYREF
#define __BOTHKNOWNFAMILYREF

#include <vector>
#include "VCF.h"

class VCFFamily;
class VCFBothKnown;
class KnownFamily;
class OptionSmall;


//////////////////// BothKnownFamilyRef ////////////////////

namespace BothKnownFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGeno *phased_vcf,
					const std::vector<std::vector<int>>& ref_haps,
					const std::vector<const KnownFamily *>& families,
					const OptionSmall& op);
};
#endif
