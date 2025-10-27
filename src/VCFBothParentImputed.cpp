#include "../include/VCFBothParentImputed.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFBothParentImputed::VCFBothParentImputed(
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const Map& map_, double w, const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf),
										records(rs),
										prog_imputer(records, map_, w) { }

VCFBothParentImputed::~VCFBothParentImputed() {
	Common::delete_all(records);
}

void VCFBothParentImputed::impute() {
	for(size_t ic = 0; ic < num_progenies(); ++ic) {
		prog_imputer.impute(ic);
	}
}
