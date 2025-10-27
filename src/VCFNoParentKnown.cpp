#include <cmath>
#include "../include/VCFNoParentKnown.h"
#include "../include/common.h"

using namespace std;

VCFNoParentKnown::VCFNoParentKnown(const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w,
									const VCFSmall *vcf) :
				VCFFamilyBase(s, vcf), records(rs),
				mat_imputer(records, true, ref_hs, map_, w),
				pat_imputer(records, false, ref_hs, map_, w),
				prog_imputer(records, map_, w) { }

VCFNoParentKnown::~VCFNoParentKnown() {
	Common::delete_all(records);
}

void VCFNoParentKnown::impute() {
	mat_imputer.impute();
	pat_imputer.impute();
	for(size_t j = 0; j < num_samples() - 2; ++j) {
		prog_imputer.impute(j);
	}
}
