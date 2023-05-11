#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFFillable.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFHeteroHomoPP ////////////////////

VCFHeteroHomoPP::VCFHeteroHomoPP(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFFillableRecord *> rs, const Map& m) :
			VCFFamily(h, s, vector<VCFFamilyRecord *>(rs.begin(), rs.end())),
			fi_records(rs), genetic_map(m) { }

bool VCFHeteroHomoPP::is_mat_hetero() const {
	return fi_records.front()->is_hetero(0);
}

string VCFHeteroHomoPP::make_seq(size_t i) const {
	stringstream	ss;
	for(auto p = this->fi_records.begin(); p != this->fi_records.end(); ++p) {
		VCFFillableRecord	*record = *p;
		const string	mat_GT = record->get_GT(0);
		const string	pat_GT = record->get_GT(1);
		const int	prog_gt = record->get_int_gt(i+2);
		if(prog_gt == -1) {
			ss << 'N';
			continue;
		}
		
		const int	mat1 = mat_GT.c_str()[0] - '0';
//		const int	mat2 = mat_GT.c_str()[2] - '0';
		const int	pat1 = pat_GT.c_str()[0] - '0';
		const int	pat2 = pat_GT.c_str()[2] - '0';
		if(pat1 == pat2) {	// mat hetero
			const int	diff = prog_gt - pat1;
			if(diff == -1 || diff == 2)
				ss << 'N';
			else
				ss << (mat1 == diff ? '0' : '1');
		}
		else {				// pat hetero
			const int	diff = prog_gt - mat1;
			if(diff == -1 || diff == 2)
				ss << 'N';
			else
				ss << (pat1 == diff ? '0' : '1');
		}
	}
	return ss.str();
}

string VCFHeteroHomoPP::impute_sample_seq(size_t i,
								const vector<double>& cMs, double min_c) {
	const string	seq = this->make_seq(i);
	if(is_all_same_without_N(seq))
		return create_same_color_string(seq);
	
	const vector<char>	hidden_states = { '0', '1' };
	const vector<char>	states = Imputer::create_states(seq);
	const string	hidden_seq = Imputer::impute(seq,
												hidden_states, states, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

double VCFHeteroHomoPP::cM(size_t i) const {
	return genetic_map.bp_to_cM(fi_records[i]->pos());
}

bool VCFHeteroHomoPP::is_all_same_without_N(const string& seq) {
	char	c = '.';
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			if(c == '.')	// initial
				c = *p;
			else if(*p != c)
				return false;
		}
	}
	return true;
}

string VCFHeteroHomoPP::create_same_color_string(const string& seq) {
	char	c = '0';	// dummy
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			c = *p;
			break;
		}
	}
	return std::string(seq.size(), c);
}

string VCFHeteroHomoPP::update_each(size_t i, size_t j, char c) {
	VCFFillableRecord	*record = this->fi_records[i];
	const int	k = (c - '0') * 2;
	if(record->is_mat_hetero())
		return record->get_GT(0).substr(k, 1) + record->get_GT(1).substr(1, 2);
	else
		return record->get_GT(0).substr(0, 2) + record->get_GT(1).substr(k, 1);
}

void VCFHeteroHomoPP::update(size_t i, const STRVEC& seqs) {
	for(size_t j = 2; j < this->samples.size(); ++j) {
		const char	c = seqs[j-2].c_str()[i];
		this->fi_records[i]->set_GT(j, this->update_each(i, j, c));
	}
}

void VCFHeteroHomoPP::impute() {
	if(this->size() == 0)
		return;
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->cM(i));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 0; i < this->num_progenies(); ++i) {
		imputed_seqs.push_back(impute_sample_seq(i, cMs, 1.0));
	}
	
	for(size_t k = 0; k < this->size(); ++k) {
		this->update(k, imputed_seqs);
	}
}

pair<ParentComb, FillType> VCFHeteroHomoPP::classify_record(
												VCFFamilyRecord *record) {
	const int	i = record->is_mat_hetero() ? 0 : 1;
	const int	j = record->is_pat_hetero() ? 0 : 1;
	if(i == 0 && j == 0) {
		return make_pair(ParentComb::P01x01, FillType::IMPUTABLE);
	}
	else if(i == 1 && j == 0) {
		if(record->mat_gt().c_str()[0] == '0')
			return make_pair(ParentComb::P00x01, FillType::PAT);
		else
			return make_pair(ParentComb::P01x11, FillType::PAT);
	}
	else if(i == 0 && j == 1) {
		if(record->pat_gt().c_str()[0] == '0')
			return make_pair(ParentComb::P00x01, FillType::MAT);
		else
			return make_pair(ParentComb::P01x11, FillType::MAT);
	}
	else {
		const char	mat_gt = record->mat_gt().c_str()[0];
		const char	pat_gt = record->pat_gt().c_str()[0];
		// FILLEDとするので、ここでimputeもしておく
		record->impute_homohomo();
		if(mat_gt == '0' && pat_gt == '0')
			return make_pair(ParentComb::P00x00, FillType::FILLED);
		else if(mat_gt == '1' && pat_gt == '1')
			return make_pair(ParentComb::P11x11, FillType::FILLED);
		else
			return make_pair(ParentComb::P00x11, FillType::FILLED);
	}
}

map<FillType, vector<VCFFillableRecord *>> VCFHeteroHomoPP::classify_records(
									const vector<VCFFamilyRecord *>& records) {
	// ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
	map<FillType, vector<VCFFillableRecord *>>	rss;
	for(size_t index = 0; index < records.size(); ++index) {
		VCFFamilyRecord	*record = records[index];
		const auto	pair = classify_record(record);
		VCFFillableRecord	*new_record = new VCFFillableRecord(record->get_v(),
						record->get_samples(), index, pair.second, pair.first);
		rss[pair.second].push_back(new_record);
	}
	return rss;
}

VCFFillable *VCFHeteroHomoPP::merge_vcf(const VCFHeteroHomoPP *mat_vcf,
					const VCFHeteroHomoPP *pat_vcf,
					const vector<VCFFillableRecord *>& homohomo_records,
					const vector<VCFFillableRecord *>& heterohetero_records) {
	for(auto p = homohomo_records.begin(); p != homohomo_records.end(); ++p) {
		VCFFillableRecord	*record = *p;
		const STRVEC&	v = record->get_v();
		const string	GT = v[9].substr(0, 1) + "|" + v[10].substr(0, 1);
		for(size_t c = 11; c != v.size(); ++c)
			record->set_GT(c-9, GT);
	}
	
	const vector<VCFFillableRecord *>&	mat_records = mat_vcf->get_fi_records();
	const vector<VCFFillableRecord *>&	pat_records = pat_vcf->get_fi_records();
	vector<VCFFillableRecord *>	records(mat_records.begin(), mat_records.end());
	records.insert(records.begin(), pat_records.begin(), pat_records.end());
	records.insert(records.begin(), homohomo_records.begin(),
										homohomo_records.end());
	records.insert(records.begin(), heterohetero_records.begin(),
										heterohetero_records.end());
	std::sort(records.begin(), records.end(), 
				[](const VCFFillableRecord *lh, const VCFFillableRecord *rh)
				{ return lh->pos() < rh->pos(); });
	return new VCFFillable(mat_vcf->get_header(),
							mat_vcf->get_samples(), records);
}

VCFFillable *VCFHeteroHomoPP::impute_by_parents(const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const STRVEC& samples, const Map& gmap) {
	VCFFamily	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
														orig_vcf, samples);
	// ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
	auto	rss = VCFHeteroHomoPP::classify_records(vcf->get_family_records());
	auto	*mat_vcf = new VCFHeteroHomoPP(vcf->get_header(),
											vcf->get_samples(),
											rss[FillType::MAT], gmap);
	auto	*pat_vcf = new VCFHeteroHomoPP(vcf->get_header(),
											vcf->get_samples(),
											rss[FillType::PAT], gmap);
	delete vcf;
	mat_vcf->impute();
	pat_vcf->impute();
	auto	new_vcf = VCFHeteroHomoPP::merge_vcf(mat_vcf, pat_vcf,
													rss[FillType::FILLED],
													rss[FillType::IMPUTABLE]);
	new_vcf->phase_hetero_hetero();
	
	// Recordは使いまわししているので、空にしてVCFだけ消す
	mat_vcf->clear_records();
	pat_vcf->clear_records();
	delete mat_vcf;
	delete pat_vcf;
	
	return new_vcf;
}
