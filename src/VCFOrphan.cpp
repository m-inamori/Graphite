#include <cmath>
#include "../include/VCFOrphan.h"
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

VCFOrphan::VCFOrphan(const STRVEC& s, const std::vector<GenoRecord *>& rs,
						const std::vector<std::vector<int>>& ref_hs,
						const Map& map_, double w, const VCFSmall *vcf) :
				VCFGenoBase(s, vcf),
				records(rs),
				imputer(records, ref_hs, map_, w) { }

VCFOrphan::~VCFOrphan() {
	Common::delete_all(records);
}

void VCFOrphan::impute(size_t io) {
	imputer.impute(io);
}
