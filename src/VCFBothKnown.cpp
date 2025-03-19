#include <cmath>
#include "../include/VCFBothKnown.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFBothKnown::VCFBothKnown(const std::vector<STRVEC>& header,
									const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				mat_imputer(new ParentImputer(records, true, ref_hs, map_, w)),
				pat_imputer(new ParentImputer(records, false, ref_hs, map_, w)),
				prog_imputer(new ProgenyImputer(records, map_, w)) { }

VCFBothKnown::~VCFBothKnown() {
	Common::delete_all(records);
	delete mat_imputer;
	delete pat_imputer;
	delete prog_imputer;
}

void VCFBothKnown::impute() {
	mat_imputer->impute();
	pat_imputer->impute();
	for(size_t j = 0; j < num_samples() - 2; ++j) {
		prog_imputer->impute(j);
	}
}
