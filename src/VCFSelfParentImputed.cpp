#include "../include/VCFSelfParentImputed.h"
#include "../include/common.h"

using namespace std;

VCFSelfParentImputed::VCFSelfParentImputed(const STRVEC& s,
											const vector<GenoRecord *>& rs,
											const Map& map_, double w,
											const VCFSmall *vcf) :
									VCFGenoBase(s, vcf),
									records(rs),
									prog_imputer(records, map_, w)
									{ }

VCFSelfParentImputed::~VCFSelfParentImputed() {
	Common::delete_all(records);
}

void VCFSelfParentImputed::impute() {
	for(size_t iprog = 0; iprog != num_progenies(); ++iprog) {
		prog_imputer.impute(iprog);
	}
}
