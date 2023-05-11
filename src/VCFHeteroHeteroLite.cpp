#include "../include/common.h"
#include "../include/VCFHeteroHeteroLite.h"

using namespace std;


//////////////////// VCFHeteroHeteroLiteRecord ////////////////////


//////////////////// VCFHeteroHeteroLite ////////////////////

VCFHeteroHeteroLite::VCFHeteroHeteroLite(const vector<STRVEC>& h,
					const STRVEC& s, vector<VCFHeteroHeteroLiteRecord *> rs) :
			VCFFamily(h, s, vector<VCFFamilyRecord *>(rs.begin(), rs.end())),
			hehe_records(rs) { }
