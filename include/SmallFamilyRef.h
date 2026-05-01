#ifndef __SMALLFAMILYREF
#define __SMALLFAMILYREF

#include <vector>
#include <string>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class KnownFamily;
class SampleManager;
class OptionSmall;


//////////////////// SmallFamilyRef ////////////////////

namespace SmallFamilyRef {
	VCFGeno *merge_vcf(VCFGeno *imputed_vcf, const VCFGenoBase *vcf,
									const std::vector<std::string>& samples);
	VCFGeno *impute_vcf_by_both_imputed_parents(const VCFSmall *orig_vcf,
												VCFGeno *phased_vcf,
												VCFGeno *imputed_vcf,
												SampleManager *sample_man,
												const OptionSmall& op_small);
	VCFGeno *impute_vcf_by_imputed_and_known_parent(
								const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	VCFGeno *impute_vcf_by_both_known_parents(
								const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	VCFGeno *impute_self_vcf(
						const VCFSmall *orig_vcf, const VCFGeno *phased_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						VCFGeno *imputed_vcf,
						const OptionSmall& op_small,
						SampleManager *sample_man);
	VCFGeno *impute_vcf_by_imputed_parent(const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	VCFGeno *impute_vcf_by_known_parent(
								const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	VCFGeno *impute_self_non_imputed_vcf(
								const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	
	std::vector<std::vector<std::string>> collect_imputed_progenies(
									std::vector<const KnownFamily *>& families,
									SampleManager *sample_man);
	// Impute families whose progenies have been imputed
	VCFGeno *impute_vcf_by_progenies(
								const VCFSmall *orig_vcf,
								const VCFGeno *phased_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								VCFGeno *imputed_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man);
	
	VCFGeno *impute_orphan_samples(
						const VCFSmall *orig_vcf,
						const std::vector<std::vector<int>>& ref_haps,
						VCFGeno *phased_vcf,
						const OptionSmall& op_small,
						SampleManager *sample_man);
	
	VCFGeno *impute_non_imputed_samples(const VCFSmall *orig_vcf,
										VCFGeno *merged_vcf,
										const OptionSmall& op_small,
										SampleManager *sample_man);
	
	// Imputes families by parent/child status until no more can be imputed
	// Updates merged_vcf in each loop with newly imputed family data
	// Returns final VCFGeno with all imputed samples included
	VCFGeno *impute(const VCFSmall *orig_vcf, VCFGeno *merged_vcf,
						const VCFGeno *ref_vcf,
						const OptionSmall& op_small,
						SampleManager *sample_man,
						bool imputes_isolated_samples);
}
#endif
