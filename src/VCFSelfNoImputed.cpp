#include <cmath>
#include "../include/VCFSelfNoImputed.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfNoImputed ////////////////////

VCFSelfNoImputed::VCFSelfNoImputed(const STRVEC& s,
									const std::vector<GenoRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w,
									const VCFSmall *vcf) :
										VCFSelfImputable(s, vcf),
										records(rs),
										imputer(records, ref_hs, map_, w) { }

VCFSelfNoImputed::~VCFSelfNoImputed() {
	Common::delete_all(records);
}

void VCFSelfNoImputed::impute() {
	imputer.impute();
}
