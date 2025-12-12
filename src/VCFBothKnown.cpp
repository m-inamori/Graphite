#include <cmath>
#include "../include/VCFBothKnown.h"
#include "../include/common.h"

using namespace std;

VCFBothKnown::VCFBothKnown(const STRVEC& s,
						   const std::vector<VCFFamilyRecord *>& rs,
						   const std::vector<std::vector<int>>& ref_haps1,
						   const std::vector<std::vector<int>>& ref_haps2,
						   const Map& map_, double w, const VCFSmall *vcf) :
								VCFFamilyBase(s, vcf), records(rs),
								mat_imputer(records, true, ref_haps1, map_, w),
								pat_imputer(records, false, ref_haps2, map_, w),
								prog_imputer(records, map_, w) { }

VCFBothKnown::~VCFBothKnown() {
	Common::delete_all(records);
}

void VCFBothKnown::impute() {
	mat_imputer.impute();
	pat_imputer.impute();
	for(size_t j = 0; j < num_samples() - 2; ++j) {
		prog_imputer.impute(j);
	}
}
