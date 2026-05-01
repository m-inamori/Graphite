#include <cstddef>

#include "../include/VCFLargeNoRef.h"

using namespace std;


//////////////////// VCFLargeNoRef ////////////////////

void VCFLargeNoRef::impute() {
	mat_imputer.impute();
	pat_imputer.impute();
	for(size_t i = 0; i < num_progenies(); ++i)
		prog_imputer.impute(i);
}
