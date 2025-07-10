#include <sstream>
#include "../include/VCFImpHeteroHomo.h"
#include "../include/Map.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFImpHeteroHomo ////////////////////

string VCFImpHeteroHomo::update_each(size_t i, size_t j, char c) const {
	
	const auto&	v = get_record(i)->get_v();
	const size_t	k = (size_t)(c - '0') * 2;
	stringstream	ss;
	if(is_mat_hetero)
		ss << v[9].c_str()[k] << v[10].substr(1);
	else
		ss << v[9].substr(0, 2) << v[10].c_str()[k];
	return ss.str();
}

void VCFImpHeteroHomo::update(size_t i, const vector<string>& seqs) {
	auto	record = get_record(i);
	const string	gt = record->get_gt(non_imputed_index());
	stringstream	ss;
	ss << gt.c_str()[0] << '|' << gt.c_str()[2];
	record->set_GT(non_imputed_index(), ss.str());
	for(size_t j = 2; j < samples.size(); ++j)
		record->set_GT(j, update_each(i, j, seqs[j-2].c_str()[i]));
}

char VCFImpHeteroHomo::determine_haplotype(int which_zero,
							int homo_int_gt, int prog_int_gt) const {
	if(prog_int_gt == -1)
		return 'N';
	
	const int	pat_int = homo_int_gt / 2;
	for(int i = 0; i < 2; ++i) {
		const int	mat_int = (i + which_zero) & 1;
		if(mat_int + pat_int == prog_int_gt)
			return '0' + i;
		else if(mat_int == 1 && mat_int + pat_int < prog_int_gt)
			return '0' + i;
		else if(mat_int == 0 && mat_int + pat_int > prog_int_gt)
			return '0' + i;
	}
	return 'N';
}

string VCFImpHeteroHomo::make_seq(size_t i) const {
	stringstream	ss;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const auto	*record = *p;
		const int	gt = record->get_int_gt(i + 2);
		if(is_mat_hetero) {
			const int	which_zero = (int)(record->mat_gt().c_str()[0] - '0');
			const int	pat_int_gt = record->pat_int_gt();
			ss << determine_haplotype(which_zero, pat_int_gt, gt);
		}
		else {
			const int	which_zero = (int)(record->pat_gt().c_str()[0] - '0');
			const int	mat_int_gt = record->mat_int_gt();
			ss << determine_haplotype(which_zero, mat_int_gt, gt);
		}
	}
	return ss.str();
}

string VCFImpHeteroHomo::impute_sample_seq(size_t j,
								const vector<double>& cMs,
								double min_c) const {
	const string	seq = make_seq(j);
	if(Imputer::is_all_same_without_N(seq))
		return Imputer::create_same_color_string(seq, '0');
	
	const string	hidden_seq = Imputer::impute(seq, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

void VCFImpHeteroHomo::impute() {
	if(records.empty())
		return;
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->cM(records[i]->pos()));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 0; i < this->num_progenies(); ++i) {
		imputed_seqs.push_back(impute_sample_seq(i, cMs, 1.0));
	}
	
	for(size_t i = 0; i < size(); ++i)
		update(i, imputed_seqs);
}
