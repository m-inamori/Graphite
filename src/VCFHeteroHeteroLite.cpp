#include "../include/common.h"
#include "../include/VCFHeteroHeteroLite.h"

using namespace std;


//////////////////// VCFHeteroHeteroLiteRecord ////////////////////


//////////////////// VCFHeteroHeteroLite ////////////////////

VCFHeteroHeteroLite::VCFHeteroHeteroLite(
							const vector<STRVEC>& h, const STRVEC& s,
							const vector<VCFHeteroHeteroLiteRecord *>& rs) :
											VCFBase(h, s), VCFSmallBase(),
											VCFFamilyBase(), records(rs) { }
