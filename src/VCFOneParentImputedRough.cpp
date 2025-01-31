#include <cmath>
#include "../include/VCFOneParentImputedRough.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputedRough::VCFOneParentImputedRough(
							const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat, const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				parent_imputer(new ParentImputer(rs, !is_mat, ref_hs, map_, w)),
				prog_imputer(new ProgenyImputer(rs, map_, w)) { }

VCFOneParentImputedRough::~VCFOneParentImputedRough() {
	Common::delete_all(records);
	delete parent_imputer;
	delete prog_imputer;
}

void VCFOneParentImputedRough::impute() {
	parent_imputer->impute();
	for(size_t ic = 0; ic < num_samples() - 2; ++ic) {
		prog_imputer->impute(ic);
	}
}
