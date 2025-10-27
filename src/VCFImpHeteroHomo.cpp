#include <sstream>
#include "../include/VCFImpHeteroHomo.h"
#include "../include/Map.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFImpHeteroHomo ////////////////////

int VCFImpHeteroHomo::update_each(size_t i, size_t j, char c) const {
	const GenoRecord	*record = this->get_record(i);
	const size_t	k = (size_t)(c - '0');
	if(is_mat_hetero) {
		const int	a1 = record->get_allele(0, k);
		const int	a2 = record->unphased(1) / 2;
		return Genotype::from_alleles(a1, a2);
	}
	else {
		const int	a1 = record->unphased(0) / 2;
		const int	a2 = record->get_allele(1, k);
		return Genotype::from_alleles(a1, a2);
	}
}

void VCFImpHeteroHomo::update(size_t i, const vector<string>& seqs) {
	auto	record = this->get_record(i);
	const int	a = record->unphased(non_imputed_index()) / 2;
	record->set_geno(non_imputed_index(), Genotype::from_alleles(a, a));
	for(size_t j = 2; j < num_samples(); ++j)
		record->set_geno(j, update_each(i, j, seqs[j-2].c_str()[i]));
}

char VCFImpHeteroHomo::determine_haplotype(int which_zero,
							int homo_int_gt, int prog_int_gt) const {
	if(Genotype::is_NA(prog_int_gt))
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
		const int	gt = record->unphased(i + 2);
		if(is_mat_hetero) {
			const int	which_zero = record->get_mat_allele(0);
			const int	pat_int_gt = record->unphased_pat();
			ss << determine_haplotype(which_zero, pat_int_gt, gt);
		}
		else {
			const int	which_zero = record->get_pat_allele(0);
			const int	mat_int_gt = record->unphased_mat();
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
		cMs.push_back(this->cM(records[i]->get_pos()));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 0; i < this->num_progenies(); ++i) {
		imputed_seqs.push_back(impute_sample_seq(i, cMs, 1.0));
	}
	
	for(size_t i = 0; i < size(); ++i)
		update(i, imputed_seqs);
}
