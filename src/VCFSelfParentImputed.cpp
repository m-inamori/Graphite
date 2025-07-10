#include "../include/SelfProgenyImputer.h"
#include "../include/VCFSelfParentImputed.h"
#include "../include/common.h"

using namespace std;

VCFSelfParentImputed::VCFSelfParentImputed(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const vector<VCFRecord *>& rs,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				imputers(create_imputers(map_, w))
				{ }

vector<SelfProgenyImputer *> VCFSelfParentImputed::create_imputers(
													const Map& map_, double w) {
	vector<SelfProgenyImputer *>	imputers(num_progenies());
	for(size_t i = 0; i < num_progenies(); ++i) {
		imputers[i] = new SelfProgenyImputer(records, i, map_, w);
	}
	return imputers;
}

VCFSelfParentImputed::~VCFSelfParentImputed() {
	Common::delete_all(records);
	Common::delete_all(imputers);
}

void VCFSelfParentImputed::impute() {
	for(auto p = imputers.begin(); p != imputers.end(); ++p) {
		(*p)->impute();
	}
}
