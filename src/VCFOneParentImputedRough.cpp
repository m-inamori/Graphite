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
				VCFOneParentImputedBase(header, s, rs), ref_records(rs),
				parent_imputer(new ParentImputer(
										ref_records, !is_mat, ref_hs, map_, w)),
				prog_imputer(new ProgenyImputer(ref_records, map_, w)) { }

VCFOneParentImputedRough::~VCFOneParentImputedRough() {
	Common::delete_all(ref_records);
	delete parent_imputer;
	delete prog_imputer;
}

void VCFOneParentImputedRough::impute() {
	parent_imputer->impute();
	for(size_t ic = 0; ic < num_samples() - 2; ++ic) {
		prog_imputer->impute(ic);
	}
}

size_t VCFOneParentImputedRough::amount() const {
	const size_t	N = num_progenies();
	const size_t	M = ref_records.size();
	const size_t	NH = parent_imputer->num_ref_haps();
	const size_t	R = (NH*NH << (N*2)) * (2*NH + 2*N - 1);
	return R * M;
}
