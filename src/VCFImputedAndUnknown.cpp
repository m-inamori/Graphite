#include <cmath>
#include "../include/VCFImputedAndUnknown.h"
#include "../include/common.h"

using namespace std;

VCFImputedAndUnknown::VCFImputedAndUnknown(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat_imputed,
							const Map& map_, double w, const VCFSmall *vcf) :
					VCFImputable(s, vcf), records(rs), ref_haps(ref_hs),
					prog_imputer(records, ref_hs, is_mat_imputed, map_, w),
					parent_imputer(records, ref_hs, is_mat_imputed, map_, w),
					progs_imputer(records, map_, w) { }

VCFImputedAndUnknown::~VCFImputedAndUnknown() {
	Common::delete_all(records);
}

size_t VCFImputedAndUnknown::amount() const {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH*NH*20;
	return R * M;
}

void VCFImputedAndUnknown::impute() {
	prog_imputer.impute(0);
	parent_imputer.impute();
	for(size_t j = 1; j < num_progenies(); ++j) {
		progs_imputer.impute(j);
	}
}
