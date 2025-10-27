#include <cmath>
#include "../include/VCFOneParentKnown.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFOneParentKnown::VCFOneParentKnown(
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const vector<vector<int>>& ref_haps,
							bool mat_known, const Map& map_, double w,
							const VCFSmall *vcf) :
				VCFFamilyBase(s, vcf),
				records(rs),
				is_mat_known(mat_known),
				parent1_imputer(records, is_mat_known, ref_haps, map_, w),
				parent2_imputer(records, !is_mat_known, ref_haps, map_, w),
				prog_imputer(records, map_, w) { }

VCFOneParentKnown::~VCFOneParentKnown() {
	Common::delete_all(records);
}

const string& VCFOneParentKnown::get_known_parent() const {
	if(is_mat_known)
		return mat();
	else
		return pat();
}

void VCFOneParentKnown::impute_known_parent() {
	parent1_imputer.impute();
}

void VCFOneParentKnown::impute() {
	parent1_imputer.impute();
	parent2_imputer.impute();
	for(size_t ic = 0; ic < num_samples() - 2; ++ic) {
		prog_imputer.impute(ic);
	}
}
