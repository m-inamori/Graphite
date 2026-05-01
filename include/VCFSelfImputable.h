#ifndef __VCFSELFIMPUTABLE
#define __VCFSELFIMPUTABLE

#include "VCFFamily.h"


//////////////////// VCFSelfImputable ////////////////////

class VCFSelfImputable : public VCFGenoBase {
public:
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFSelfImputable *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFSelfImputable *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
public:
	VCFSelfImputable(const STRVEC& s, const VCFSmall *vcf) :
										VCFGenoBase(s, vcf) { }
	virtual ~VCFSelfImputable() { }
	
	std::size_t num_progenies() const { return num_samples() - 1; }
	
	// amount of computation
	virtual std::size_t amount() const = 0;
	virtual void impute() = 0;
	
public:
	static void impute_in_thread(void *config);
	static void impute_VCFs(std::vector<VCFSelfImputable *>& vcfs, int T);
};
#endif
