#include <cmath>
#include "../include/VCFSelfNoImputed.h"
#include "../include/common.h"

using namespace std;

VCFSelfNoImputed::VCFSelfNoImputed(const STRVEC& s,
									const std::vector<GenoRecord *>& rs,
									const std::vector<std::vector<int>>& ref_hs,
									const Map& map_, double w,
									const VCFSmall *vcf) :
										VCFGenoBase(s, vcf),
										VCFMeasurable(map_),
										records(rs),
										imputer(records, ref_hs, map_, w) { }

VCFSelfNoImputed::~VCFSelfNoImputed() {
	Common::delete_all(records);
}

void VCFSelfNoImputed::impute() {
	imputer.impute();
}
