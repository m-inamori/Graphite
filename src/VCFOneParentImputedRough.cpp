#include <cmath>
#include "../include/VCFOneParentImputedRough.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputedRough::VCFOneParentImputedRough(const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									bool is_mat, const Map& map_, double w,
									const VCFSmall *vcf) :
						VCFOneParentImputedBase(s, rs, vcf),
						records(rs),
						parent_imputer(records, !is_mat, ref_hs, map_, w),
						prog_imputer(records, map_, w) { }

VCFOneParentImputedRough::~VCFOneParentImputedRough() {
	Common::delete_all(records);
}

void VCFOneParentImputedRough::impute() {
	parent_imputer.impute();
	for(size_t ic = 0; ic < num_samples() - 2; ++ic) {
		prog_imputer.impute(ic);
	}
}

size_t VCFOneParentImputedRough::amount() const {
	const size_t	N = num_progenies();
	const size_t	M = records.size();
	const size_t	NH = parent_imputer.num_ref_haps();
	const size_t	R = (NH*NH << (N*2)) * (2*NH + 2*N - 1);
	return R * M;
}
