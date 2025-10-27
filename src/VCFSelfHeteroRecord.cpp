#include <sstream>
#include <stack>
#include "../include/VCFSelfHetero.h"
#include "../include/common.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"
#include "../include/inverse_graph.h"
#include "../include/option.h"

using namespace std;


//////////////////// VCFSelfHeteroRecord ////////////////////

vector<int> VCFSelfHeteroRecord::progeny_gts() const {
	const vector<int>&	geno = this->get_geno();
	return vector<int>(geno.begin() + 1, geno.end());
}

void VCFSelfHeteroRecord::set_haplo(int h) {
	this->set_geno(0, h == 0 ? Genotype::PH_01 : Genotype::PH_10);
}

void VCFSelfHeteroRecord::set_int_gt_by_which_comes_from(int w1, int w2,
																	size_t i) {
	const int	parent_gt = this->get_geno()[0];
	
	const int	gt = Genotype::from_alleles((parent_gt >> w1) & 1,
											(parent_gt >> w2) & 1);
	set_geno(i + 1, gt);
}
