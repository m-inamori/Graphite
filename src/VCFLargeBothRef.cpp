#include <cstddef>

#include "../include/VCFLargeBothRef.h"

using namespace std;


//////////////////// VCFLargeBothRef ////////////////////

void VCFLargeBothRef::impute() {
	for(size_t i = 0; i < num_progenies(); ++i)
		imputer.impute(i);
}
