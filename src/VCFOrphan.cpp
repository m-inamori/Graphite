#include <cmath>
#include "../include/VCFOrphan.h"
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

VCFOrphan::VCFOrphan(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				imputer(new OrphanImputer(records, ref_hs, map_, w)) { }

VCFOrphan::~VCFOrphan() {
	Common::delete_all(records);
	delete imputer;
}

void VCFOrphan::impute(size_t io) {
	imputer->impute(io);
}
