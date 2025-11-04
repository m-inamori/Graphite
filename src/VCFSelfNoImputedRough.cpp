#include "../include/VCFSelfNoImputedRough.h"
#include "../include/common.h"

using namespace std;

VCFSelfNoImputedRough::VCFSelfNoImputedRough(const STRVEC& s,
											const vector<GenoRecord *>& rs,
											const vector<vector<int>>& ref_hs,
											const Map& map_, double w,
											const VCFSmall *vcf) :
									VCFGenoBase(s, vcf),
									VCFMeasurable(map_),
									records(rs),
									parent_imputer(records, ref_hs, map_, w),
									prog_imputer(records, map_, w)
									{ }

VCFSelfNoImputedRough::~VCFSelfNoImputedRough() {
	Common::delete_all(records);
}

void VCFSelfNoImputedRough::impute() {
	parent_imputer.impute();
	for(size_t iprog = 0; iprog != num_progenies(); ++iprog) {
		prog_imputer.impute(iprog);
	}
}
