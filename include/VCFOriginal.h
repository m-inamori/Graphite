#ifndef __VCFORIGINAL
#define __VCFORIGINAL

#include "VCF.h"

class PedigreeTable;


//////////////////// VCFOriginal ////////////////////

class VCFOriginal : public VCFHuge {
public:
	VCFOriginal(const std::vector<STRVEC>& h, const STRVEC& s, VCFReader *r) :
															VCFHuge(h, s, r) { }
	VCFOriginal(const VCFOriginal&) = delete;
	VCFOriginal& operator=(const VCFOriginal&) = delete;
	~VCFOriginal() { }
	
public:
	static VCFOriginal *read(const std::string& path);
};
#endif
