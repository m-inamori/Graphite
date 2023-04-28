#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "common.h"
#include "VCFOneParentPhased.h"
#include "VCFFillable.h"
#include "Map.h"
#include "Pedigree.h"
#include "Imputer.h"

using namespace std;


//////////////////// VCFOneParentPhased ////////////////////

VCFOneParentPhased::VCFOneParentPhased(const vector<STRVEC>& h, const STRVEC& s,
					vector<VCFFamilyRecord *> rs, bool mat_p, const Map& m) :
			VCFFamily(h, s, vector<VCFFamilyRecord *>(rs.begin(), rs.end())),
			is_mat_phased(mat_p), genetic_map(m) { }

bool VCFOneParentPhased::is_mat_hetero() const {
	return family_records.front()->is_hetero(0);
}

char VCFOneParentPhased::determine_which_comes_from(
									VCFFamilyRecord *record, size_t i) const {
	if(record->is_NA(0) || record->is_NA(1) || record->is_NA(i))
		return 'N';
		
	const string	mat_GT = record->get_GT(0);
	const string	pat_GT = record->get_GT(1);
	const int		prog_gt = record->get_int_gt(i);
	const int		mat1 = mat_GT.c_str()[0] - '0';
	const int		mat2 = mat_GT.c_str()[2] - '0';
	const int		pat1 = pat_GT.c_str()[0] - '0';
	const int		pat2 = pat_GT.c_str()[2] - '0';
	if(is_mat_phased) {
		if(!(mat1 != mat2 && pat1 == pat2))
			return 'N';
		
		const int	diff = prog_gt - pat1;
		if(diff == -1 || diff == 2)
			return 'N';
		else
			return diff == mat1 ? '0' : '1';
	}
	else {	// pat is phased
		if(!(mat1 == mat2 && pat1 != pat2))
			return 'N';
		
		const int	diff = prog_gt - mat1;
		if(diff == -1 || diff == 2)
			return 'N';
		else
			return diff == pat1 ? '0' : '1';
	}
}

string VCFOneParentPhased::make_seq(size_t i) const {
	stringstream	ss;
	for(auto p = family_records.begin(); p != family_records.end(); ++p) {
		VCFFamilyRecord	*record = *p;
		ss << this->determine_which_comes_from(record, i);
	}
	return ss.str();
}

double VCFOneParentPhased::cM(size_t i) const {
	return genetic_map.bp_to_cM(records[i]->pos());
}

string VCFOneParentPhased::impute_sample_seq(size_t i,
								const vector<double>& cMs, double min_c) const {
	const string	seq = make_seq(i);
	return Imputer::impute_seq(seq, cMs, min_c);
}

void VCFOneParentPhased::update_each(size_t i, size_t k, char c) {
	VCFRecord	*record = this->records[k];
	const STRVEC&	v = record->get_v();
	const size_t	j = (size_t)(c - '0') * 2;
	if(is_mat_phased) {
		const char	mat_gt = v[9].c_str()[j];
		const char	mat_int_gt = (int)(mat_gt - '0');
		const int	remain = record->get_int_gt(i) - mat_int_gt;
		if(remain <= 0) {
			record->set_GT(i, string(1, mat_gt) + "|0");
			record->set_GT(1, record->get_int_gt(1) == 0 ? "0|0" : "0|1");
		}
		else {
			record->set_GT(i, string(1, mat_gt) + "|1");
			record->set_GT(1, record->get_int_gt(1) == 2 ? "1|1" : "1|0");
		}
	}
	else {
		const char	pat_gt = v[10].c_str()[j];
		const int	pat_int_gt = (int)(pat_gt - '0');
		const int	remain = record->get_int_gt(i) - pat_int_gt;
		if(remain <= 0) {
			record->set_GT(i, string("0|") + pat_gt);
			record->set_GT(0, record->get_int_gt(0) == 0 ? "0|0" : "0|1");
		}
		else {
			record->set_GT(i, string("1|") + pat_gt);
			record->set_GT(0, record->get_int_gt(0) == 2 ? "1|1" : "1|0");
		}
	}
}

void VCFOneParentPhased::update(size_t i, const string& seq) {
	for(size_t k = 0; k < this->size(); ++k) {
		this->update_each(i, k, seq.c_str()[k]);
	}
}

void VCFOneParentPhased::impute() {
	if(this->size() == 0)
		return;
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->cM(i));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 2; i < this->num_samples(); ++i) {
		const string	imputed_seq = impute_sample_seq(i, cMs, 1.0);
		this->update(i, imputed_seq);
	}
}

VCFFamily *VCFOneParentPhased::impute_by_parent(const VCFSmall *orig_vcf,
									const VCFSmall *parent_imputed_vcf,
									const STRVEC& samples,
									bool is_mat_phased, const Map& gmap) {
	VCFFamily	*vcf = VCFFamily::create_by_two_vcfs(parent_imputed_vcf,
														orig_vcf, samples);
	vector<VCFFamilyRecord *>	records = vcf->get_family_records();
	auto	new_vcf = new VCFOneParentPhased(vcf->get_header(),
							vcf->get_samples(), records, is_mat_phased, gmap);
	new_vcf->impute();
//	vector<VCFFamilyRecord *>	new_records = new_vcf->get_family_records();
	VCFFamily	*imputed_vcf = new VCFFamily(vcf->get_header(),
//									vcf->get_samples(), new_records);
									vcf->get_samples(), records);
	new_vcf->clear_records();	// As imputed_vcf uses the records as they are
	delete new_vcf;
	vcf->clear_records();
	delete vcf;
	return imputed_vcf;
}

// samplesに対応するRecordを作る
VCFRecord *VCFOneParentPhased::merge_records(vector<VCFFamily *>& vcfs,
											size_t i, const STRVEC& samples) {
	// sampleがどのVCFのどの位置にあるか
	map<string, pair<const VCFFamily *, size_t>>	dic_pos;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const VCFFamily	*vcf = *p;
		for(size_t k = 0; k < vcf->num_samples(); ++k)
			dic_pos[vcf->get_samples()[k]] = make_pair(vcf, k + 9);
	}
	
	// 上の辞書を元にRecordを作る
	const STRVEC&	v1 = vcfs.front()->get_records()[i]->get_v();
	STRVEC	v(v1.begin(), v1.begin() + 9);
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		auto	q = dic_pos[*p];
		v.push_back(q.first->get_records()[i]->get_v()[q.second]);
	}
	return new VCFRecord(v, samples);
}
