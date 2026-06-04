#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/ParentProgenyImputer.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFOneParentImputed ////////////////////

VCFOneParentImputed::VCFOneParentImputed(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat_imp, const Map& map_,
							double w, const VCFSmall *vcf) :
					VCFImputable(s, vcf),
					records(rs),
					is_mat_imputed(is_mat_imp),
					ref_haps(ref_hs),
					imputer(rs, is_mat_imp, ref_hs, map_, w) { }

VCFOneParentImputed::~VCFOneParentImputed() {
	Common::delete_all(records);
}

size_t VCFOneParentImputed::amount() const {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH*NH * (2*NH - 1);
	return R * M;
}

STRVEC VCFOneParentImputed::imputed_samples() const {
	STRVEC	ss = samples;
	const size_t	index = is_mat_imputed ? 0 : 1;
	ss.erase(ss.begin() + index);
	return ss;
}

void VCFOneParentImputed::impute() {
	imputer.impute();
}
