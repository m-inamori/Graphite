#include <cmath>
#include "../include/VCFSelfRef.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfRef ////////////////////

VCFSelfRef::VCFSelfRef(const STRVEC& s, const std::vector<GenoRecord *>& rs,
							const Map& map_, double w, const VCFSmall *vcf) :
										VCFSelfImputable(s, vcf),
										records(rs),
										imputer(records, map_, w) { }

VCFSelfRef::~VCFSelfRef() {
	Common::delete_all(records);
}

void VCFSelfRef::impute() {
	for(size_t i = 0; i < num_progenies(); ++i) {
		imputer.impute(i);
	}
}
