#include <cstddef>

#include "../include/VCFLargeBothRef.h"

using namespace std;


//////////////////// VCFLargeBothRef ////////////////////

VCFLargeBothRef::~VCFLargeBothRef() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

STRVEC VCFLargeBothRef::imputed_samples() const {
	STRVEC	ss = samples;
	ss.erase(ss.begin(), ss.begin() + 2);
	return ss;
}

void VCFLargeBothRef::impute() {
	for(size_t i = 0; i < num_progenies(); ++i)
		imputer.impute(i);
}
