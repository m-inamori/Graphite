#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/ParentProgenyImputer.h"
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
					imputer(records, ref_hs, should_impute_mat_, map_, w),
					should_impute_mat(should_impute_mat_) { }

VCFOneParentImputed::~VCFOneParentImputed() {
	Common::delete_all(records);
}

size_t VCFOneParentImputed::amount() const {
	return ref_haps.size() << (samples.size()*2-3);
}

STRVEC VCFOneParentImputed::imputed_samples() const {
	STRVEC	ss(samples.size()-2);
	const size_t	index = should_impute_mat ? 0 : 1;
	ss[0] = samples[index];
	std::copy(samples.begin() + 3, samples.end(), ss.begin() + 1);
	return ss;
}

void VCFOneParentImputed::impute() {
	imputer.impute();
}
