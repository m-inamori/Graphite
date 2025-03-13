#include <cmath>
#include "../include/VCFProgenyImputed.h"
#include "../include/ParentImputerByProgeny.h"
#include "../include/common.h"

using namespace std;

VCFProgenyImputed::VCFProgenyImputed(const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFRecord *>& rs,
							const std::vector<std::vector<int>>& ref_haps,
							bool is_mat_known, const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(),
				records(rs),
				imputer(new ParentImputerByProgeny(records, ref_haps,
													is_mat_known, map_, w)) { }

VCFProgenyImputed::~VCFProgenyImputed() {
	delete imputer;
	Common::delete_all(records);
}

void VCFProgenyImputed::impute() {
	imputer->impute();
}
