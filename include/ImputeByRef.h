// ImputeByRef.h
#ifndef __IMPUTEBYREF
#define __IMPUTEBYREF

#include "VCF.h"

class VCFGeno;
class SampleManager;
class Map;
class Option;
class Materials;


//////////////////// ImputeByRef ////////////////////

namespace ImputeByRef {
	void display_chromosome_info(const VCFSmall *orig_vcf);
	VCFGeno *remove_reference_samples(const VCFGeno *vcf,
								const std::vector<std::string>& ref_samples);
	
	VCFGeno *impute_vcf_chr(const VCFSmall *orig_vcf, const VCFSmall *ref_vcf_,
										SampleManager *sample_man,
										const Map& gmap, const Option& Option);
	
	void impute(VCFHuge *vcf, const Materials *materials,
											const Option& option);
}
#endif
