#include "../include/VCFSelfHeteroRecord.h"

using namespace std;


//////////////////// VCFSelfHeteroRecord ////////////////////

vector<int> VCFSelfHeteroRecord::progeny_gts() const {
	const vector<int>&	geno = this->get_geno();
	return vector<int>(geno.begin() + 1, geno.end());
}

void VCFSelfHeteroRecord::set_haplo(int h) {
	this->set_geno(0, h == 0 ? Genotype::PH_01 : Genotype::PH_10);
}

