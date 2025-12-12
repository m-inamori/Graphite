#ifndef __REFERENCEHAPLOTYPE
#define __REFERENCEHAPLOTYPE

#include <vector>

class VCFGenoBase;
class VCFGeno;
class SampleManager;


namespace ReferenceHaplotype {
	std::vector<std::vector<int>> extract_haplotypes(
										const VCFGeno *phased_vcf,
										const SampleManager *sample_man);
	int count_same_alleles(const std::vector<int>& gts,
											const std::vector<int>& ref_hap);
	std::vector<std::vector<int>> filter_haplotypes(
								const std::vector<std::vector<int>>& ref_haps,
								const std::vector<int>& gts, std::size_t upper);
}
#endif
