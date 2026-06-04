#ifndef __VCFIMPUTABLE
#define __VCFIMPUTABLE

#include "VCFFamily.h"


//////////////////// VCFImputable ////////////////////

class VCFImputable : public VCFFamilyBase {
public:
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFImputable *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFImputable *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
public:
	VCFImputable(const STRVEC& s, const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf) { }
	virtual ~VCFImputable() { }
	
	// amount of computation
	virtual std::size_t amount() const = 0;
	virtual STRVEC imputed_samples() const = 0;
	virtual void impute() = 0;
	
public:
	static void impute_in_thread(void *config);
	static void impute_VCFs(std::vector<VCFImputable *>& vcfs, int T);
	// Merge imputed sample columns from multiple VCFs into a single VCF
	// in the specified sample order.
	static VCFGeno *join(const std::vector<VCFImputable *>& vcfs,
											const STRVEC& samples);
};
#endif
