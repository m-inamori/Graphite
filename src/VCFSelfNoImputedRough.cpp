#include "../include/SelfParentImputerLessImputed.h"
#include "../include/SelfProgenyImputer.h"
#include "../include/VCFSelfNoImputedRough.h"
#include "../include/common.h"

using namespace std;

VCFSelfNoImputedRough::VCFSelfNoImputedRough(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				imputer(new SelfParentImputerLessImputed(records,
														ref_hs, map_, w)),
				imputers(create_imputers(map_, w))
				{ }

vector<SelfProgenyImputer *> VCFSelfNoImputedRough::create_imputers(
													const Map& map_, double w) {
	vector<SelfProgenyImputer *>	imputers(num_progenies());
	for(size_t i = 0; i < num_progenies(); ++i) {
		imputers[i] = new SelfProgenyImputer(records, i, map_, w);
	}
	return imputers;
}

VCFSelfNoImputedRough::~VCFSelfNoImputedRough() {
	Common::delete_all(records);
	delete imputer;
	Common::delete_all(imputers);
}

void VCFSelfNoImputedRough::impute() {
	imputer->impute();
	for(auto p = imputers.begin(); p != imputers.end(); ++p) {
		(*p)->impute();
	}
}
