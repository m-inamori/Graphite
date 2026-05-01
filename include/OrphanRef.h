#ifndef __ORPHANREF
#define __ORPHANREF

#include "VCF.h"

class VCFGenoBase;
class VCFGeno;
class OptionSmall;


//////////////////// OrphanRef ////////////////////

namespace OrphanRef {
	VCFGenoBase *impute(const std::vector<std::string>& samples,
							const VCFSmall *orig_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const VCFGeno *phased_vcf,
							const OptionSmall& op);
};
#endif
