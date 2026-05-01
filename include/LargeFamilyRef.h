#ifndef __LARGEFAMILYREF
#define __LARGEFAMILYREF

#include <vector>

class VCFSmall;
class VCFGeno;
class VCFFamilyBase;
class VCFImputable;
class KnownFamily;
class Map;
class Option;


//////////////////// LargeFamilyRef ////////////////////

namespace LargeFamilyRef {
	VCFImputable *create_family_vcf(const VCFFamilyBase *vcf_family,
									const VCFGeno *ref_vcf, const Map& gmap);
	VCFGeno *impute(const std::vector<const KnownFamily *>& families,
					const VCFSmall *orig_vcf,
					const VCFGeno *ref_vcf,
					const Map& gmap, const Option& option);
}
#endif
