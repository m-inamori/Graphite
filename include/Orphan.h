#ifndef __ORPHAN
#define __ORPHAN

#include <vector>
#include "VCF.h"

class VCFGenoBase;
class VCFOrphan;
class VCFOrphanRough;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// Orphan ////////////////////

namespace Orphan {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		VCFOrphan	*vcf;
		
		ConfigThread(std::size_t i, std::size_t n, VCFOrphan *vcf_) :
									first(i), num_threads(n), vcf(vcf_) { }
	};
	
	void impute_small_in_thread(void *config);
	void impute_small_VCF(VCFOrphan *vcf, int T);
	
	struct ConfigThreadRough {
		std::size_t	first;
		std::size_t	num_threads;
		VCFOrphanRough	*vcf;
		
		ConfigThreadRough(std::size_t i, std::size_t n, VCFOrphanRough *vcf_) :
										first(i), num_threads(n), vcf(vcf_) { }
	};
	
	void impute_small_in_thread_rough(void *config);
	void impute_small_VCF_rough(VCFOrphanRough *vcf, int T);
	
	VCFGenoBase *impute(const std::vector<std::string>& samples,
							const VCFSmall *orig_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const OptionSmall& op);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
													const OptionSmall& op);
	// upper NH which passes is_small
	std::size_t compute_upper_NH(std::size_t M, const OptionSmall& op);
};
#endif
