#include <cstddef>

#include "../include/VCFLargeOneRef.h"

using namespace std;


//////////////////// VCFLargeOneRef ////////////////////

VCFLargeOneRef::~VCFLargeOneRef() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

STRVEC VCFLargeOneRef::imputed_samples() const {
	STRVEC	ss = samples;
	const size_t	index = is_mat_ref ? 0 : 1;
	ss.erase(ss.begin() + index);
	return ss;
}

void VCFLargeOneRef::impute() {
	parent_imputer.impute();
	for(size_t i = 0; i < num_progenies(); ++i)
		prog_imputer.impute(i);
}
