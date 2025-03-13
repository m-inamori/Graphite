#ifndef __VCFORIGINAL
#define __VCFORIGINAL

#include "VCF.h"

class PedigreeTable;


//////////////////// VCFOriginal ////////////////////

class VCFOriginal : public VCFHuge {
public:
	VCFOriginal(const std::vector<STRVEC>& h, const STRVEC& s, VCFReader *r) :
															VCFHuge(h, s, r) { }
	
public:
	static VCFOriginal *read(const std::string& path);
};
#endif
