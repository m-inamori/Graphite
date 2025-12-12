#ifndef __VCFIMPUTABLE
#define __VCFIMPUTABLE

#include "VCFFamily.h"


//////////////////// VCFImputable ////////////////////

class VCFImputable : public VCFFamilyBase {
public:
	VCFImputable(const STRVEC& s, const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf) { }
	virtual ~VCFImputable() { }
	
	// amount of computation
	virtual std::size_t amount() const = 0;
	virtual void impute() = 0;
};
#endif
