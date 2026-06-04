#include <cmath>
#include "../include/VCFParentsKnownProgenyImputed.h"
#include "../include/ParentImputer.h"
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

VCFParentsKnownProgenyImputed::VCFParentsKnownProgenyImputed(
							const STRVEC& s,
							const vector<VCFFamilyRecord *>& rs,
							const vector<vector<int>>& ref_haps_mat,
							const vector<vector<int>>& ref_haps_pat,
							const Map& map_, double w, const VCFSmall *vcf) :
				VCFImputable(s, vcf),
				records(rs),
				NH(ref_haps_mat.size()),
				imputer(records, ref_haps_mat, ref_haps_pat, map_, w) { }

VCFParentsKnownProgenyImputed::~VCFParentsKnownProgenyImputed() {
	Common::delete_all(records);
}

void VCFParentsKnownProgenyImputed::impute() {
	imputer.impute();
}

size_t VCFParentsKnownProgenyImputed::amount() const {
	return NH * NH * (8 * NH + 4);
}
