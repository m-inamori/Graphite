#ifndef __SMALLFAMILY
#define __SMALLFAMILY

#include <vector>
#include "../include/VCF.h"

class VCFSmall;
class VCFFamily;
class VCFFillable;
class SampleManager;
class KnownFamily;
class Map;
class Option;

namespace SmallFamily {
	VCFSmall *impute_vcf_by_both_imputed_parents(
						const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
						SampleManager *sample_man,
						const Map& geno_map, const Option *option);
	
	VCFSmall *impute_vcf_by_imputed_and_known_parent(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const Map& geno_map,
							SampleManager *sample_man, const Option *option);
	
	VCFSmall *impute_vcf_by_both_known_parents(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const Map& gmap,
							SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_vcf_by_imputed_parent(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	VCFSmall *impute_one_parent_vcf(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	std::vector<std::vector<std::string>> collect_imputed_progenies(
									std::vector<const KnownFamily *>& families,
									SampleManager *sample_man);
	// Impute families whose progenies have been imputed
	VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_vcf_by_known_parent(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_orphan_samples(const VCFSmall *orig_vcf,
								VCFSmall *merged_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);

	
	std::vector<std::vector<int>> extract_haplotypes(
										const VCFSmall *phased_vcf,
										const SampleManager *sample_man);
	
	VCFSmall *impute_small_family_VCFs(const VCFSmall *orig_vcf,
										VCFSmall *merged_vcf,
										const Map& geno_map,
										SampleManager *sample_man,
										const Option *option);
	
	VCFSmall *impute_vcf_by_progenies_core(const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const Map& geno_map,
							SampleManager *sample_man, int num_threads);
	VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
	
	VCFSmall *impute_iolated_samples(
					const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
					SampleManager *sample_man, const STRVEC& samples,
					const Map& gmap, bool modify_genotypes, int num_threads);
}
#endif
