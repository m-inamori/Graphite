#ifndef __IMPUTENOREF
#define __IMPUTENOREF

class VCFSmall;
class VCFHuge;
class VCFGeno;
class Materials;
class SampleManager;
class Map;
class Option;

namespace ImputeNoRef {
	void display_chromosome_info(const VCFSmall *orig_vcf);
	VCFGeno *impute_vcf_chr(const VCFSmall *orig_vcf, SampleManager *sample_man,
									const Map& geno_map, const Option& option);
	void impute(VCFHuge *vcf, const Materials *materials, const Option& option);
}
#endif
