#include <sstream>
#include <cassert>

#include "../include/VCFSelfFillableRecord.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfFillableRecord ////////////////////

VCFSelfFillableRecord *VCFSelfFillableRecord::copy() const {
	return new VCFSelfFillableRecord(v, samples, index, type, comb, probs);
}

string VCFSelfFillableRecord::gt_from_parent(int mat_from, int pat_from) const {
	stringstream	ss;
	ss << v[9].c_str()[mat_from*2-2] << '|' << v[9].c_str()[pat_from*2-2];
	return ss.str();
}

int VCFSelfFillableRecord::find_geno_type(const string& type) const {
	const auto	vec = Common::split(this->v[8], ':');
	for(auto p = vec.begin(); p != vec.end(); ++p) {
		if(*p == type)
			return p - vec.begin();
	}
	return -1;
}

int VCFSelfFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	int j = is_mat ? 0 : 1;
	const string&	parent_gt = this->get_gt(0);
	return parent_gt.c_str()[0] == this->get_gt(i).c_str()[j*2] ? 1 : 2;
}

VCFSelfFillableRecord *VCFSelfFillableRecord::convert(
									const VCFImpSelfRecord *record) {
	const SelfFillType	type = record->get_fill_type();
	const auto&	probs = record->parse_PL();
	return new VCFSelfFillableRecord(record->get_v(), record->get_samples(),
						 record->get_index(), type, record->get_comb(), probs);
}

// after imputed
int VCFSelfFillableRecord::from_which_chrom(const VCFSelfFillableRecord *record,
													size_t i, bool is_mat) {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}
