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
				VCFImputable(s, vcf),
				records(rs),
				is_mat_known(mat_known),
				parent_imputer(records, is_mat_known, ref_haps, map_, w) { }

VCFOneParentKnown::~VCFOneParentKnown() {
	Common::delete_all(records);
}

void VCFOneParentKnown::impute() {
	parent_imputer.impute();
}

STRVEC VCFOneParentKnown::imputed_samples() const {
	const size_t	index = phased_index();
	return { samples[index] };
}
