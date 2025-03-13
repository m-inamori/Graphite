#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputedBase::VCFOneParentImputedBase(
							const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs) :
				VCFBase(header, s), VCFSmallBase(), records(rs) { }
