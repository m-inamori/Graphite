#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFOneParentImputed ////////////////////

VCFOneParentImputed::VCFOneParentImputed(const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									bool should_impute_mat_, const Map& map_,
									double w, const VCFSmall *vcf) :
						VCFImputable(s, vcf),
						records(rs),
						ref_haps(ref_hs),
						imputer(records, should_impute_mat_, ref_hs, map_, w),
						should_impute_mat(should_impute_mat_) { }

VCFOneParentImputed::~VCFOneParentImputed() {
	Common::delete_all(records);
}

size_t VCFOneParentImputed::amount() const {
	const size_t	NH = ref_haps.size();
	const size_t	N = samples.size() - 2;
	return NH*NH * (NH+2*N-1) << (2*N);
}

STRVEC VCFOneParentImputed::imputed_samples() const {
	STRVEC	ss(samples.size()-1);
	const size_t	index = should_impute_mat ? 0 : 1;
	ss[0] = samples[index];
	std::copy(samples.begin() + 2, samples.end(), ss.begin() + 1);
	return ss;
}

void VCFOneParentImputed::impute() {
	imputer.impute();
}
