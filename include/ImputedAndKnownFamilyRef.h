#ifndef __ONEPHASEDFAMILYREF
#define __ONEPHASEDFAMILYREF

#include <array>
#include "VCFImpFamilyRecord.h"
#include "VCFFillableRecord.h"
#include "TypeDeterminer.h"

class VCFOneParentImputed;
class VCFSmallFillable;
class VCFImputable;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ImputedAndKnownFamilyRef ////////////////////

namespace ImputedAndKnownFamilyRef {
	VCFGeno *impute(const VCFSmall *orig_vcf,
						const VCFGeno *phased_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						const std::vector<const KnownFamily *>& families,
						const STRVEC& non_imputed_parents,
						const OptionSmall& op);
};
#endif
