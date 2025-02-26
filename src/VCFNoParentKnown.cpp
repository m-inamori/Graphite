#include <cmath>
#include "../include/VCFNoParentKnown.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFNoParentKnown::VCFNoParentKnown(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									bool is_mat_known,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				parent_imputer1(new ParentImputer(records, is_mat_known,
															ref_hs, map_, w)),
				parent_imputer2(new ParentImputer(records, !is_mat_known,
															ref_hs, map_, w)),
				prog_imputer(new ProgenyImputer(records, map_, w)) { }

VCFNoParentKnown::~VCFNoParentKnown() {
	Common::delete_all(records);
	delete parent_imputer1;
	delete parent_imputer2;
	delete prog_imputer;
}

void VCFNoParentKnown::impute() {
	parent_imputer1->impute();
	parent_imputer2->impute();
	for(size_t j = 0; j < num_samples() - 2; ++j) {
		prog_imputer->impute(j);
	}
}
