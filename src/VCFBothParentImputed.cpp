#include <cmath>
#include "../include/VCFBothParentImputed.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFBothParentImputed::VCFBothParentImputed(
							const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), records(rs),
				prog_imputer(new ProgenyImputer(rs, map_, w)) { }

VCFBothParentImputed::~VCFBothParentImputed() {
	Common::delete_all(records);
	delete prog_imputer;
}

void VCFBothParentImputed::impute() {
	for(size_t ic = 0; ic < num_samples() - 2; ++ic) {
		prog_imputer->impute(ic);
	}
}
