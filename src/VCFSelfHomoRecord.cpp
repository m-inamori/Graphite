#include "../include/common.h"
#include "../include/VCFSelfHomoRecord.h"

using namespace std;


//////////////////// VCFSelfHomoRecord ////////////////////

vector<int> VCFSelfHomoRecord::gts() const {
	if(this->comb == ParentComb::P00x00)
		return vector<int>(num_samples(), Genotype::PH_00);
	else
		return vector<int>(num_samples(), Genotype::PH_11);
}

VCFSelfHomoRecord *VCFSelfHomoRecord::impute() const {
	return new VCFSelfHomoRecord(this->pos, this->gts(), this->index,
											this->wrong_type, this->comb);
}
