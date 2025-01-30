#include <cmath>
#include "../include/VCFNoParentImputed.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFNoParentImputed::VCFNoParentImputed(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w) :
					VCFBase(header, s), VCFSmallBase(), records(rs),
					mat_imputer(new ParentImputer(rs, true, ref_hs, map_, w)),
					pat_imputer(new ParentImputer(rs, false, ref_hs, map_, w)),
					prog_imputer(new ProgenyImputer(rs, map_, w)) { }

VCFNoParentImputed::~VCFNoParentImputed() {
	Common::delete_all(records);
	delete mat_imputer;
	delete pat_imputer;
	delete prog_imputer;
}

void VCFNoParentImputed::impute() {
	mat_imputer->impute();
	pat_imputer->impute();
	for(size_t j = 0; j < num_samples() - 2; ++j) {
		prog_imputer->impute(j);
	}
}
