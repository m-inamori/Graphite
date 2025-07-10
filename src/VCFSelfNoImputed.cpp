#include <cmath>
#include "../include/VCFSelfNoImputed.h"
#include "../include/SelfImputer.h"
#include "../include/common.h"

using namespace std;

VCFSelfNoImputed::VCFSelfNoImputed(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				imputer(new SelfImputer(records, ref_hs, map_, w)) { }

VCFSelfNoImputed::~VCFSelfNoImputed() {
	Common::delete_all(records);
	delete imputer;
}

void VCFSelfNoImputed::impute() {
	imputer->impute();
}
