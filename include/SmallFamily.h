#ifndef __SMALLFAMILY
#define __SMALLFAMILY

#include <vector>
#include "../include/VCF.h"

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFFamily;
class VCFFillable;
class SampleManager;
class KnownFamily;
class Map;
class Option;
class OptionSmall;

namespace SmallFamily {
	VCFGeno *impute_vcf_by_both_imputed_parents(
						const VCFSmall *orig_vcf, const VCFGenoBase *merged_vcf,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_vcf_by_imputed_and_known_parent(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_vcf_by_both_known_parents(
						const VCFSmall *orig_vcf, const VCFGeno *imputed_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_vcf_by_imputed_parent(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	VCFGeno *impute_self_vcf(const VCFSmall *orig_vcf,
								const VCFGeno *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								SampleManager *sample_man,
								const OptionSmall& op_small);
	VCFGeno *impute_self_non_imputed_vcf(const VCFSmall *orig_vcf,
								const VCFGeno *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								SampleManager *sample_man,
								const OptionSmall& op_small);
	VCFGeno *impute_one_parent_vcf(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	std::vector<std::vector<std::string>> collect_imputed_progenies(
									std::vector<const KnownFamily *>& families,
									SampleManager *sample_man);
	// Impute families whose progenies have been imputed
	VCFGeno *impute_vcf_by_progenies(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_vcf_by_known_parent(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_orphan_samples(
						const VCFSmall *orig_vcf, VCFGeno *merged_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	std::vector<std::vector<int>> extract_haplotypes(
										const VCFGeno *phased_vcf,
										const SampleManager *sample_man);
	
	VCFGeno *impute_vcf_by_progenies(
						const VCFSmall *orig_vcf, const VCFGeno *merged_vcf,
						SampleManager *sample_man, const OptionSmall& op_small);
	
	VCFGeno *impute_small_family(
						const VCFSmall *orig_vcf, VCFGeno *merged_vcf,
						const Map& geno_map, const Option *option,
						SampleManager *sample_man);
	
	VCFGeno *impute_non_imputed_samples(
						const VCFSmall *orig_vcf, VCFGeno *merged_vcf,
						SampleManager *sample_man, const OptionSmall& op_small);
}
#endif
