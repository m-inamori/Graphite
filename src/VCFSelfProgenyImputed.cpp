#include "../include/SelfParentImputer.h"
#include "../include/SelfProgenyImputer.h"
#include "../include/VCFSelfProgenyImputed.h"
#include "../include/common.h"

using namespace std;

VCFSelfProgenyImputed::VCFSelfProgenyImputed(const STRVEC& s,
											 const vector<GenoRecord *>& rs,
											 const vector<vector<int>>& ref_hs,
											 size_t ic_, const Map& map_,
											 double w, const VCFSmall *vcf) :
				VCFGenoBase(s, vcf), records(rs), ic(ic_),
				imputer(records, ref_hs, ic, map_, w),
				imputers(create_imputers(map_, w))
				{ }

vector<SelfProgenyImputer *> VCFSelfProgenyImputed::create_imputers(
													const Map& map_, double w) {
	vector<SelfProgenyImputer *>	imputers(num_progenies());
	for(size_t i = 0; i < num_progenies(); ++i) {
		imputers[i] = new SelfProgenyImputer(records, i, map_, w);
	}
	return imputers;
}

VCFSelfProgenyImputed::~VCFSelfProgenyImputed() {
	Common::delete_all(records);
	Common::delete_all(imputers);
}

void VCFSelfProgenyImputed::impute() {
	imputer.impute();
	for(auto p = imputers.begin(); p != imputers.end(); ++p) {
		(*p)->impute();
	}
}
