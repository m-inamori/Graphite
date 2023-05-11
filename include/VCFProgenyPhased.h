#ifndef __VCFPROGENYPHASED
#define __VCFPROGENYPHASED

#include "VCFFamily.h"


//////////////////// VCFOneParentPhased ////////////////////

class VCFProgenyPhased : public VCFFamily {
	const std::vector<std::size_t>	phased_progeny_indices;
	
public:
	VCFProgenyPhased(const std::vector<STRVEC>& h, const STRVEC& s,
										std::vector<VCFFamilyRecord *> rs,
										const std::vector<std::size_t>& ppi);
	~VCFProgenyPhased() { }
	
	void determine_parent(bool is_mat);
	void impute();
	
public:
	static VCFProgenyPhased *impute_by_progeny(const VCFSmall *orig_vcf,
										const VCFSmall *imputed_vcf,
										const STRVEC& samples,
										const std::vector<std::size_t>& ppi);
};
#endif
