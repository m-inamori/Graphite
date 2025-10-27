#include <sstream>
#include <cassert>

#include "../include/VCFSelfFillableRecord.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfFillableRecord ////////////////////

int VCFSelfFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	int j = is_mat ? 0 : 1;
	const int	parent_gt = this->get_geno()[0];
	return (parent_gt & 1) == ((this->get_geno()[i] >> j) & 1) ? 1 : 2;
}

VCFSelfFillableRecord *VCFSelfFillableRecord::convert(
									const VCFImpSelfRecord *record,
									const vector<VCFRecord::Probs>& probs) {
	const SelfFillType	type = record->get_fill_type();
	return new VCFSelfFillableRecord(record->get_pos(), record->get_geno(),
						 record->get_index(), type, record->get_comb(), probs);
}

// after imputed
int VCFSelfFillableRecord::from_which_chrom(const VCFSelfFillableRecord *record,
													size_t i, bool is_mat) {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}
