#ifndef __BOTHIMPUTEDFAMILYREF
#define __BOTHIMPUTEDFAMILYREF

#include <vector>
#include "VCF.h"

class VCFFamily;
class VCFBothParentImputed;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// BothImputedFamilyRef ////////////////////

namespace BothImputedFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGenoBase *phased_vcf,
					const std::vector<const KnownFamily *>& families,
					const OptionSmall& op);
};
#endif
