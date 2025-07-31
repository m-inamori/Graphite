#include "../include/common.h"
#include "../include/VCFSelfHomoRecord.h"

using namespace std;


//////////////////// VCFSelfHomoRecord ////////////////////

vector<string> VCFSelfHomoRecord::gts() const {
	if(this->comb == ParentComb::P00x00)
		return vector<string>(samples.size(), "0|0");
	else
		return vector<string>(samples.size(), "1|1");
}

VCFSelfHomoRecord *VCFSelfHomoRecord::impute() const {
	vector<string>	w = v;
	const vector<string>	gts1 = this->gts();
	for(size_t i = 9; i < w.size(); ++i)
		w[i] = gts1[i-9] + w[i].substr(3);
	return new VCFSelfHomoRecord(w, this->samples, this->index,
									this->wrong_type, this->comb);
}
