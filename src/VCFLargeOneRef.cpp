#include <cstddef>

#include "../include/VCFLargeOneRef.h"

using namespace std;


//////////////////// VCFLargeOneRef ////////////////////

VCFLargeOneRef::~VCFLargeOneRef() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFLargeOneRef::impute() {
	parent_imputer.impute();
	for(size_t i = 0; i < num_progenies(); ++i)
		prog_imputer.impute(i);
}
