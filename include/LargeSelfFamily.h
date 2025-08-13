#ifndef __LARGESELFFAMILY
#define __LARGESELFFAMILY

#include <vector>
#include "VCF.h"
#include "ClassifyRecord.h"

class VCFSmall;
class VCFSelfFillable;
class VCFSelfHeteroRecord;
class VCFImpSelfRecord;
class Family;
class KnownFamily;
class SampleManager;
class Map;
class Option;

namespace LargeSelfFamily {
	std::pair<std::vector<VCFSelfHeteroRecord *>,
			  std::vector<VCFImpSelfRecord *>>
						divide_records(const VCFSmall *vcf, const Option *op);
	VCFSmall *extract_parents(const std::vector<VCFSelfFillable *>& vcfs);
	VCFSmall *impute(const VCFSmall *orig_vcf, VCFSmall *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map, const Option *op);
}
#endif
