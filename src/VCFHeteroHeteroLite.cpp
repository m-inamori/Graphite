#include "../include/common.h"
#include "../include/VCFHeteroHeteroLite.h"

using namespace std;


//////////////////// VCFHeteroHeteroLiteRecord ////////////////////


//////////////////// VCFHeteroHeteroLite ////////////////////

VCFHeteroHeteroLite::VCFHeteroHeteroLite(const STRVEC& s,
							const vector<VCFHeteroHeteroLiteRecord *>& rs,
							const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf),
										records(rs) { }
