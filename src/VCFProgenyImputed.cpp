#include <cmath>
#include "../include/VCFProgenyImputed.h"
#include "../include/common.h"

using namespace std;

VCFProgenyImputed::VCFProgenyImputed(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_haps,
							bool is_mat_known, const Map& map_,
							double w, const VCFSmall *vcf) :
				VCFFamilyBase(s, vcf),
				records(rs),
				imputer(records, ref_haps, is_mat_known, map_, w) { }

VCFProgenyImputed::~VCFProgenyImputed() {
	Common::delete_all(records);
}

void VCFProgenyImputed::impute() {
	imputer.impute();
}
