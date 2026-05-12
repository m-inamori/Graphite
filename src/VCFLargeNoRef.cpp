#include <cstddef>

#include "../include/VCFLargeNoRef.h"

using namespace std;


//////////////////// VCFLargeNoRef ////////////////////

VCFLargeNoRef::~VCFLargeNoRef() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFLargeNoRef::impute() {
	mat_imputer.impute();
	pat_imputer.impute();
	for(size_t i = 0; i < num_progenies(); ++i)
		prog_imputer.impute(i);
}
