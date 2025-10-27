#ifndef __LARGESELFFAMILY
#define __LARGESELFFAMILY

#include <vector>
#include "VCF.h"
#include "ClassifyRecord.h"

class VCFSmall;
class VCFGeno;
class VCFSelfFillable;
class VCFSelfHeteroRecord;
class VCFImpSelfRecord;
class KnownFamily;
class Map;
class Option;

namespace LargeSelfFamily {
	std::pair<std::vector<VCFSelfHeteroRecord *>,
			  std::vector<VCFImpSelfRecord *>>
						divide_records(const VCFGeno *vcf, const Option *op);
	VCFGeno *extract_parents(const std::vector<VCFSelfFillable *>& vcfs);
	VCFGeno *impute(const VCFSmall *orig_vcf, VCFGeno *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map, const Option *op);
}
#endif
