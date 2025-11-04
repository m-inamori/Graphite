#include "../include/VCFSelfProgenyImputed.h"
#include "../include/common.h"

using namespace std;

VCFSelfProgenyImputed::VCFSelfProgenyImputed(const STRVEC& s,
											 const vector<GenoRecord *>& rs,
											 const vector<vector<int>>& ref_hs,
											 size_t ic_, const Map& map_,
											 double w, const VCFSmall *vcf) :
								VCFGenoBase(s, vcf), records(rs), ic(ic_),
								parent_imputer(records, ref_hs, ic, map_, w),
								prog_imputer(records, map_, w)
								{ }

VCFSelfProgenyImputed::~VCFSelfProgenyImputed() {
	Common::delete_all(records);
}

void VCFSelfProgenyImputed::impute() {
	parent_imputer.impute();
	for(size_t iprog = 0; iprog != num_progenies(); ++iprog) {
		if(iprog != ic)
			prog_imputer.impute(iprog);
	}
}
