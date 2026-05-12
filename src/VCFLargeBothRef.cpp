#include <cstddef>

#include "../include/VCFLargeBothRef.h"

using namespace std;


//////////////////// VCFLargeBothRef ////////////////////

VCFLargeBothRef::~VCFLargeBothRef() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFLargeBothRef::impute() {
	for(size_t i = 0; i < num_progenies(); ++i)
		imputer.impute(i);
}
