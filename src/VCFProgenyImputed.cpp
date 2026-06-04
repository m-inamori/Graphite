#include <cmath>
#include "../include/VCFProgenyImputed.h"
#include "../include/common.h"

using namespace std;

VCFProgenyImputed::VCFProgenyImputed(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_haps,
							bool should_impute_mat_, const Map& map_,
							double w, const VCFSmall *vcf) :
				VCFImputable(s, vcf),
				records(rs),
				should_impute_mat(should_impute_mat_),
				imputer(records, ref_haps, should_impute_mat, map_, w) { }

VCFProgenyImputed::~VCFProgenyImputed() {
	Common::delete_all(records);
}

STRVEC VCFProgenyImputed::imputed_samples() const {
	if(should_impute_mat)
		return { samples[0] };
	else
		return { samples[1] };
}

void VCFProgenyImputed::impute() {
	imputer.impute();
}
