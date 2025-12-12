#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/ParentProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputed::VCFOneParentImputed(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat_imputed, const Map& map_,
							double w, const VCFSmall *vcf) :
					VCFImputable(s, vcf),
					records(rs),
					ref_haps(ref_hs),
					imputer(rs, is_mat_imputed, ref_hs, map_, w) { }

VCFOneParentImputed::~VCFOneParentImputed() {
	Common::delete_all(records);
}

size_t VCFOneParentImputed::amount() const {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH*NH * (2*NH - 1);
	return R * M;
}

void VCFOneParentImputed::impute() {
	imputer.impute();
}
