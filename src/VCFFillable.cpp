#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>
#include <cassert>
#include "../include/VCFFillable.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// Genotype ////////////////////

Genotype::Genotype(const string& s) : gt1(s.c_str()[0]), gt2(s.c_str()[2]),
												phasing(s.c_str()[1] == '|') { }

pair<char,char> Genotype::gts() const {
	return pair<char,char>(gt1, gt2);
}

bool Genotype::includes(char gt) const {
	return gt == gt1 || gt == gt2;
}

bool Genotype::conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
												bool considers_phasing) {
	if(considers_phasing && this->phasing) {
		return mat_gt.includes(gt1) && pat_gt.includes(gt2);
	}
	else {
		if(gt1 == gt2)
			return mat_gt.includes(gt1) && pat_gt.includes(gt2);
		else
			return (mat_gt.includes(gt1) && pat_gt.includes(gt2)) ||
					(pat_gt.includes(gt1) && mat_gt.includes(gt2));
	}
}

bool Genotype::is_valid(const string& gt, int mat_gt, int pat_gt) {
	if(gt.length() < 3)
		return false;
	
	const string	mat_gts = possible_gts(mat_gt);
	const string	pat_gts = possible_gts(pat_gt);
	return (mat_gts.find(gt.substr(0, 1)) != string::npos &&
			pat_gts.find(gt.substr(2, 1)) != string::npos) ||
		   (mat_gts.find(gt.substr(2, 1)) != string::npos &&
			pat_gts.find(gt.substr(0, 1)) != string::npos);
}

string Genotype::possible_gts(int gt) {
	switch(gt) {
		case  0: return "0";
		case  3: return "1";
		default: return "01";
	}
}

int Genotype::sum_gt(const string& gt) {
	return (int)((gt.c_str()[0] - '0') + (gt.c_str()[2] - '0'));
}

bool Genotype::is_all_NA(const vector<string>& GTs) {
	for(auto p = GTs.begin(); p != GTs.end(); ++p) {
		if(*p != "./.")
			return false;
	}
	return true;
}


//////////////////// VCFFillableRecord ////////////////////

#define P(a,b) make_pair(a,b)

vector<pair<int, int>> VCFFillableRecord::possible_phasings() const {
	switch(this->comb) {
		case ParentComb::P00x00: return { P(0, 0) };
		case ParentComb::P00x01: return { P(0, 1), P(0, 2), P(1, 0), P(2, 0),
										  P(1, 1), P(1, 2), P(2, 1), P(2, 2) };
		case ParentComb::P01x01: return { P(1, 1), P(1, 2), P(2, 1), P(2, 2),
										  P(0, 1), P(0, 2), P(1, 0), P(2, 0),
										  P(1, 3), P(2, 3), P(3, 1), P(3, 2) };
		case ParentComb::P00x11: return { P(0, 3), P(3, 0) };
		case ParentComb::P01x11: return { P(1, 3), P(2, 3), P(3, 1), P(3, 2),
										  P(1, 1), P(1, 2), P(2, 1), P(2, 2) };
		case ParentComb::P11x11: return { P(3, 3) };
		default: break;
	}
	
	// any pair is possible
	vector<pair<int, int>>	pairs;
	for(int i = 0; i < 16; ++i)
		pairs.push_back(P(i >> 2, i & 3));
	return pairs;
}

STRVEC VCFFillableRecord::prog_gts() const {
	STRVEC	vec(v.size() - 11);
	std::copy(v.begin() + 11, v.end(), vec.begin());
	return vec;
}

VCFFillableRecord *VCFFillableRecord::copy() const {
	return new VCFFillableRecord(v, samples, index, type, comb);
}

tuple<int,int,int> VCFFillableRecord::count_gts() const {
	const auto	gts = get_progeny_int_gts();
	int	counter[3] = { 0, 0, 0 };
	for(auto p = gts.begin() + 2; p != gts.end(); ++p) {
		if(*p != -1)
			counter[*p] += 1;
	}
	return make_tuple(counter[0], counter[1], counter[2]);
}

string VCFFillableRecord::gt_from_parent(int mat_from, int pat_from) const {
	stringstream	ss;
	ss << v[9].c_str()[mat_from*2-2] << '|' << v[10].c_str()[pat_from*2-2];
	return ss.str();
}

string VCFFillableRecord::gt_from_mat(int mat_from, int c) const {
	const string	gt = this->get_GT(c - 9);
	const int		int_gt = this->get_int_gt(c - 9);
	const char		mat_gt = this->v[9].c_str()[mat_from*2-2];
	const string	pat_GT = this->get_GT(1);
	if(int_gt == 0) {
		if(mat_gt == '0')
			return "0|0";
		else if(pat_GT.find('0') != string::npos)
			return "1|0";
		else
			return "1|1";
	}
	else if(int_gt == 2) {
		if(mat_gt == '1')
			return "1|1";
		else if(pat_GT.find('1') != string::npos)
			return "0|1";
		else
			return "0|0";
	}
	else if(int_gt == 1) {
		if(mat_gt == '0')
			return "0|1";
		else
			return "1|0";
	}
	else {
		return ".|.";
	}
}

string VCFFillableRecord::gt_from_pat(int pat_from, int c) const {
	const string	gt = this->get_GT(c - 9);
	const int		int_gt = this->get_int_gt(c - 9);
	const char		pat_gt = this->v[10].c_str()[pat_from*2-2];
	const string	mat_GT = this->get_GT(0);
	if(int_gt == 0) {
		if(pat_gt == '0')
			return "0|0";
		else if(mat_GT.find('0') != string::npos)
			return "0|1";
		else
			return "1|1";
	}
	else if(int_gt == 2) {
		if(pat_gt == '1')
			return "1|1";
		else if(mat_GT.find('1') != string::npos)
			return "1|0";
		else
			return "0|0";
	}
	else if(int_gt == 1) {
		if(pat_gt == '0')
			return "1|0";
		else
			return "1|0";
	}
	else {
		return ".|.";
	}
}

vector<vector<double>> VCFFillableRecord::make_probability_table() const {
	// 全ての親のGenotypeの組合せで確率を求める
	// (0/0, 0/0), (0/0, 0/1), (0/0, 1/1), (0/1, 0/1), (0/1, 1/1), (1/1, 1/1)
	const vector<double>	ps = {1.0, 0.5, 0.0};
	vector<vector<double>>	pss_;
	for(int i = 0; i < 3; ++i) {
		const double	mat = ps[i];
		for(int j = i; j < 3; ++j) {
			const double	pat = ps[j];
			const vector<double>	v = {mat*pat, mat+pat-2*mat*pat,
													(1.0-mat)*(1.0-pat)};
			pss_.push_back(v);
		}
	}
	
	vector<vector<double>>	pss(6, vector<double>(3));
	const double	p_miss = 0.01;
	for(int i = 0; i < 6; ++i) {
		for(int j = 0; j < 3; ++j)
			pss[i][j] = (pss_[i][j] + p_miss) / (1.0 + 3 * p_miss);
	}
	return pss;
}

void VCFFillableRecord::phase() {
	// 両親ともホモになっている前提
	const string	gt_mat = this->v[9].substr(0, 1);
	const string	gt_pat = this->v[10].substr(0, 1);
	this->set_mat_GT(gt_mat + "|" + gt_mat);
	this->set_pat_GT(gt_pat + "|" + gt_pat);
	for(size_t i = 2; i < this->num_samples(); ++i)
		this->set_GT(i, gt_mat + "|" + gt_pat);
}

int VCFFillableRecord::mat_from(int c) const {
	if(this->v[c].c_str()[0] == '.')
		return 0;
	else if(!this->is_hetero(0))
		return 0;
	else if(this->v[c].c_str()[0] == this->v[9].c_str()[0])
		return 1;
	else
		return 2;
}

int VCFFillableRecord::pat_from(int c) const {
	if(this->v[c].c_str()[0] == '.')
		return 0;
	else if(!this->is_hetero(1))
		return 0;
	else if(this->v[c].c_str()[2] == this->v[10].c_str()[0])
		return 1;
	else
		return 2;
}

int VCFFillableRecord::find_geno_type(const string& type) const {
	const auto	vec = Common::split(this->v[8], ':');
	for(auto p = vec.begin(); p != vec.end(); ++p) {
		if(*p == type)
			return p - vec.begin();
	}
	return -1;
}

void VCFFillableRecord::fill_PGT() {
	const int	i_GT = this->find_geno_type("GT");
	assert(i_GT != -1);
	const int	i_PGT = this->find_geno_type("PGT");
	if(i_PGT == -1)
		return;
	
	for(size_t j = 9U; j < this->v.size(); ++j) {
		auto	vec = Common::split(this->v[j], ':');
		vec[i_PGT] = v[i_GT];
		this->v[j] = Common::join(vec, ':');
	}
}

string VCFFillableRecord::inverse_gt(const string& gt, bool inv) const {
	if(inv)
		return string(gt == "0" ? "1" : "0");
	else
		return gt;
}

string VCFFillableRecord::inverse_prog_gt(const string& gt,
											bool inv_mat, bool inv_pat) const {
	return inverse_gt(gt.substr(0, 1), inv_mat) + "|" +
			inverse_gt(gt.substr(2, 1), inv_pat) + gt.substr(3);
}

STRVEC VCFFillableRecord::inverse_prog_gts(const STRVEC& prog_gts,
											bool inv_mat, bool inv_pat) const {
	STRVEC	inv_prog_gts;
	for(auto p = prog_gts.begin(); p != prog_gts.end(); ++p) {
		inv_prog_gts.push_back(inverse_prog_gt(*p, inv_mat, inv_pat));
	}
	return inv_prog_gts;
}

void VCFFillableRecord::inverse_parents_gts(bool inv_mat, bool inv_pat) {
	// both must be hetero
	if(inv_mat)
		set_GT(0, v[9].substr(2, 1) + "|" + v[9].substr(0, 1));
	if(inv_pat)
		set_GT(1, v[10].substr(2, 1) + "|" + v[10].substr(0, 1));
}

bool VCFFillableRecord::is_same_gts(const string& gt1,
									const string& gt2) const {
	if(gt2 == "0/1")
		return gt1 == "0|1" || gt1 == "1|0";
	else if(gt2 == "0/0")
		return gt1 == "0|0";
	else if(gt2 == "1/1")
		return gt1 == "1|1";
	else
		return false;
}

bool VCFFillableRecord::is_near_prog_gts(const STRVEC& gts) const {
	int	num = 0;
	int	dist = 0;
	for(size_t i = 0U; i < gts.size(); ++i) {
		if(v[i+11] != "0/1")
			num += 1;
		if(!is_same_gts(gts[i], v[i+11]))
			dist += 1;
	}
	return dist < max(1, num / 2);
}

void VCFFillableRecord::modify_gts(const STRVEC& new_prog_gts) {
	for(int i = 0; i < 4; ++i) {
		const bool	inv_mat = (i >> 1) == 1;
		const bool	inv_pat = (i & 1) == 1;
		STRVEC	inv_prog_gts = inverse_prog_gts(new_prog_gts, inv_mat, inv_pat);
		if(is_near_prog_gts(inv_prog_gts)) {
			inverse_parents_gts(inv_mat, inv_pat);
			std::copy(inv_prog_gts.begin(), inv_prog_gts.end(), v.begin() + 11);
			return;
		}
	}
	std::copy(new_prog_gts.begin(), new_prog_gts.end(), v.begin() + 11);
}

void VCFFillableRecord::modify_parents_type() {
	if(!is_00x11() &&
			((get_GT(0) == "0|0" && get_GT(1) == "1|1") ||
			 (get_GT(0) == "1|1" && get_GT(1) == "0|0")))
		comb = ParentComb::P00x11;
}

int VCFFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	int j = is_mat ? 0 : 1;
	const string&	parent_gt = this->get_gt(j);
	return parent_gt.c_str()[0] == this->get_gt(i).c_str()[j*2] ? 1 : 2;
}

void VCFFillableRecord::set(const STRVEC& new_v, FillType new_type) {
	v = new_v;
	type = new_type;
	VCFFamilyRecord::set(new_v);
}

VCFRecord *VCFFillableRecord::integrate_records(
								const vector<VCFFillableRecord *>& records) {
	for(auto p = records.begin(); p != records.end(); ++p) {
		if((*p)->is_unable())
			return NULL;
	}
	
	STRVEC	v(records.front()->v.begin(), records.front()->v.begin() + 9);
	for(auto p = records.begin(); p != records.end(); ++p) {
		v.insert(v.end(), (*p)->v.begin() + 9, (*p)->v.end());
	}
	return new VCFRecord(v, records.front()->get_samples());
}

int VCFFillableRecord::hash(int d) const {
	int	hash_ = 0;
	const string&	str_pos = this->v[1];
	for(auto p = str_pos.begin(); p != str_pos.end(); ++p) {
		hash_ = (hash_ + (*p - '0')) % d;
	}
	return hash_;
}

string VCFFillableRecord::decide_by_majority(const vector<string>& GTs) const {
	// ex) ['0/0', '0/1', '0/0', '1/1'] -> { '0/0': 2, '0/1': 1, '1/1': 1 }
	map<string, int>	counter;
	for(auto p = GTs.begin(); p != GTs.end(); ++p)
		counter[*p] += 1;
	
	// ex) { '0/0': 2, '0/1': 1, '1/1': 1 } -> { 2: ['0/0'], 1: ['0/1', '1/1'] }
	map<int, vector<string>>	dic_num;
	for(auto p = counter.begin(); p != counter.end(); ++p)
		dic_num[p->second].push_back(p->first);
	
	pair<int, vector<string>>	max_pair;
	for(auto p = dic_num.begin(); p != dic_num.end(); ++p) {
		if(p->first > max_pair.first)
			max_pair = *p;
	}
	
	vector<string>	max_GTs = max_pair.second;
	std::sort(max_GTs.begin(), max_GTs.end());	// Python版と合わせるため
	if(max_GTs.size() == 1)
		return max_GTs.front();
	else	// 最多が複数ある場合は乱数的なものを使う
		return max_GTs[this->hash((int)max_GTs.size())];
}

void VCFFillableRecord::swap_parents(int i, const string& GT) {
	if(GT == this->v[i+9].substr(0, 3))
		return;
	
	if(GT == "0|0" || GT == "1|1") {
		const bool	is_mat_00 = (i == 0) ^ (GT == "1|1");
		this->set_GT(0, is_mat_00 ? "0|0" : "1|1");
		this->set_GT(1, is_mat_00 ? "1|1" : "0|0");
		
		// 子どもも入れ換えなければならない
		const string	prog_GT = is_mat_00 ? "0|1" : "1|0";
		for(size_t i = 2; i < this->samples.size(); ++i)
			this->set_GT(i, prog_GT);
	}
}

string VCFFillableRecord::decide_duplicated_Genotype(
									const vector<VCFFillableRecord *>& records,
									const vector<pair<int, int>>& positions) {
	// 全部./.ならそのまま
	// 子どもがいたらそれ
	// ./.を除いて全部同じだったらそれ
	// ./.と0/0 x 1/1で全部なら0/0 x 1/1で多数決
	// ./.と0/0 x 1/1を除いて全部同じだったらそれ
	// ./.と0/0 x 1/1を除いて複数種あれば多数決
	STRVEC	GTs;
	for(auto p = positions.begin(); p != positions.end(); ++p)
		GTs.push_back(records[p->first]->get_GT(p->second));
	
	if(Genotype::is_all_NA(GTs))
		return "./.";
	
	// 子どもが優先
	for(size_t i = 0; i < GTs.size(); ++i) {
		const int	j = positions[i].second;
		const string&	GT = GTs[i];
		if(j >= 2 && GT != "./.")
			return GT;
	}
	
	STRVEC	GTs_less_NA;
	for(auto p = GTs.begin(); p != GTs.end(); ++p) {
		if(*p != "./.")
			GTs_less_NA.push_back(*p);
	}
	
	if(Common::is_all_same(GTs_less_NA))
		return GTs_less_NA.front();
	
	STRVEC	GTs_less_00x11;
	for(size_t i = 0; i < GTs.size(); ++i) {
		const auto	*record = records[positions[i].first];
		const string&	GT = GTs[i];
		if(GT != "./." && record->comb != ParentComb::P00x11)
			GTs_less_00x11.push_back(GT);
	}
	
	if(GTs_less_00x11.empty())	// only 0/0 x 1/1
		return records.front()->decide_by_majority(GTs_less_NA);
	else if(Common::is_all_same(GTs_less_00x11))
		return GTs_less_00x11.front();
	else
		return records.front()->decide_by_majority(GTs_less_00x11);
}

VCFFillableRecord *VCFFillableRecord::convert(
									const VCFImpFamilyRecord *record) {
	const FillType	type = record->get_fill_type();
	return new VCFFillableRecord(record->get_v(), record->get_samples(),
								 record->get_index(), type, record->get_comb());
}

VCFRecord *VCFFillableRecord::merge(const vector<VCFFillableRecord *>& records,
														const STRVEC& samples) {
	STRVEC	v = records.front()->v;
	for(auto p = records.begin() + 1; p != records.end(); ++p)
		v.insert(v.end(), (*p)->v.begin(), (*p)->v.end());
	return new VCFRecord(v, samples);
}

void VCFFillableRecord::integrate_each_sample(
									const vector<VCFFillableRecord *>& records,
									const vector<pair<int, int>>& positions) {
	const string	GT = decide_duplicated_Genotype(records, positions);
	for(auto p = positions.begin(); p != positions.end(); ++p) {
		VCFFillableRecord	*record = records[p->first];
		if(record->is_00x11() && p->second <= 1)
			record->swap_parents(p->second, GT);
	}
}

bool VCFFillableRecord::is_all_same_GT(
						const vector<VCFFillableRecord *>& records,
						const vector<pair<int, int>>& pos_samples) {
	if(pos_samples.size() == 1)
		return true;
	
	const auto	pos_sample = pos_samples.front();
	const string	GT0 = records[pos_sample.first]->get_GT(pos_sample.second);
	for(auto p = pos_samples.begin() + 1; p != pos_samples.end(); ++p) {
		const string	GT = records[p->first]->get_GT(p->second);
		if(GT != GT0)
			return false;
	}
	return true;
}

VCFRecord *VCFFillableRecord::integrate(
							const vector<VCFFillableRecord *>& records,
							const vector<string>& samples,
							const vector<vector<pair<int, int>>>& pos_samples) {
	const VCFFillableRecord	*record = records.front();
	// 0|0 x 1|1なら親のGenotypeを交換できる
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		if(!is_all_same_GT(records, *p)) {
			integrate_each_sample(records, *p);
		}
	}
	
	// 交換したあとにGenotypeを集める
	vector<string>	v(record->v.begin(), record->v.begin() + 9);
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		const pair<int, int>&	pos = p->front();
		v.push_back(records[pos.first]->get_gt(pos.second));
	}
	return new VCFRecord(v, samples);
}

// phasingされている前提
int VCFFillableRecord::from_which_chrom(const VCFFillableRecord *record,
													size_t i, bool is_mat) {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}


//////////////////// VCFFillable::RecordSet ////////////////////

string VCFFillable::RecordSet::gt_each(int i,
										const VCFFillableRecord *record) const {
	return record == NULL ? "" : record->get_gt(i);
}

string VCFFillable::RecordSet::prev_mat_gt(int i) const {
	return gt_each(i, prev_mat_record);
}

string VCFFillable::RecordSet::next_mat_gt(int i) const {
	return gt_each(i, next_mat_record);
}

string VCFFillable::RecordSet::prev_pat_gt(int i) const {
	return gt_each(i, prev_pat_record);
}

string VCFFillable::RecordSet::next_pat_gt(int i) const {
	return gt_each(i, next_pat_record);
}

int VCFFillable::RecordSet::prev_mat_from(std::size_t i) const {
	if(prev_mat_record == NULL)
		return 0;
	return prev_mat_record->from_which_chrom(i, true);
}
int VCFFillable::RecordSet::next_mat_from(std::size_t i) const {
	if(next_mat_record == NULL)
		return 0;
	return next_mat_record->from_which_chrom(i, true);
}
int VCFFillable::RecordSet::prev_pat_from(std::size_t i) const {
	if(prev_pat_record == NULL)
		return 0;
	return prev_pat_record->from_which_chrom(i, false);
}
int VCFFillable::RecordSet::next_pat_from(std::size_t i) const {
	if(next_pat_record == NULL)
		return 0;
	return next_pat_record->from_which_chrom(i, false);
}

bool VCFFillable::RecordSet::is_mat_prev_near() const {
	return record->pos() * 2 < prev_mat_record->pos() + next_mat_record->pos();
}

bool VCFFillable::RecordSet::is_pat_prev_near() const {
	return record->pos() * 2 < prev_pat_record->pos() + next_pat_record->pos();
}

int VCFFillable::RecordSet::near_mat_from(size_t i) const {
	return is_mat_prev_near() ? prev_mat_from(i) : next_mat_from(i);
}

int VCFFillable::RecordSet::near_pat_from(size_t i) const {
	return is_pat_prev_near() ? prev_pat_from(i) : next_pat_from(i);
}

VCFFillable::Pair VCFFillable::RecordSet::select_nearest_froms(
								const vector<Pair>& pairs, size_t i) const {
	if(pairs.size() == 4U) {
		return Pair(near_mat_from(i), near_pat_from(i));
	}
	else if(pairs[0].first == pairs[1].first) {			// matが同じ
		if(is_pat_prev_near())
			return Pair(pairs[0].first, prev_pat_from(i));
		else
			return Pair(pairs[0].first, next_pat_from(i));
	}
	else if(pairs[0].second == pairs[1].second) {	// patが同じ
		if(is_mat_prev_near())
			return Pair(prev_mat_from(i), pairs[0].second);
		else
			return Pair(next_mat_from(i), pairs[1].second);
	}
	else {	// 両親とも乗り換えている（滅多にない）
		return Pair(near_mat_from(i), near_pat_from(i));
	}
}

VCFFillable::Pair VCFFillable::RecordSet::select_pair(const vector<Pair>& pairs,
												size_t i, bool selected) const {
	if(pairs.empty())
		return Pair(0, 0);
	else if(pairs.size() == 1U)
		return pairs.front();
	else if(!Genotype::is_valid(record->get_gt(i), record->mat_int_gt(),
													record->pat_int_gt()))
		return select_nearest_froms(pairs, i);
	else if(selected)
		return select_nearest_froms(pairs, i);
	
	vector<Pair>	new_pairs;	// 元々のGenotypeと同じになるペア
	for(auto p = pairs.begin(); p != pairs.end(); ++p) {
		string	parent_gt = record->gt_from_parent(p->first, p->second);
		if(Genotype::sum_gt(record->get_gt(i)) == Genotype::sum_gt(parent_gt))
			new_pairs.push_back(*p);
	}
	const Pair	selected_pair = select_pair(new_pairs, i, true);
	if(selected_pair.first != 0)
		return selected_pair;
	else
		return select_pair(pairs, i, true);
}

// 親のどちらの染色体から来ているかの確率
vector<double> VCFFillable::RecordSet::probs_from_which_chrom(
									int prev_from, int next_from) const {
	static const double ps[] = {
		0.5, 0.9, 0.1, 0.9, 0.99, 0.5, 0.1, 0.5, 0.01
	};
	const double	prob1 = ps[prev_from+next_from*3];
	vector<double>	probs = { prob1, 1.0 - prob1 };
	return probs;
}

vector<double> VCFFillable::RecordSet::probs_from_which_chrom(
												size_t i, bool is_mat) const {
//	if((is_mat && record->is_mat_homo()) || (!is_mat && record->is_pat_homo()))
//		return vector<double>(2U, 0.5);
	if(is_mat)
		return probs_from_which_chrom(prev_mat_from(i), next_mat_from(i));
	else
		return probs_from_which_chrom(prev_pat_from(i), next_pat_from(i));
}

double VCFFillable::RecordSet::likelihood_each(const string& gt,
									const vector<double>& probs_mat,
									const vector<double>& probs_pat,
									int mat_phasing, int pat_phasing) const {
	const int	sum = Genotype::sum_gt(gt);
	double	likelihood = 0.0;
	for(int k = 0; k < 4; ++k) {
		const int	i = k >> 1;		// 母親は0|1か1|0か
		const int	j = k & 1;
		if((((mat_phasing >> i) & 1) + ((pat_phasing >> j) & 1)) == sum)
			likelihood += probs_mat[i] * probs_pat[j];
	}
	if(likelihood == 0.0)	// 該当する組合せが無い
		return log(0.0001);
	else
		return log(likelihood);
}

double VCFFillable::RecordSet::compute_phasing_likelihood_each(size_t i,
									int mat_phasing, int pat_phasing) const {
	const auto	probs_mat = this->probs_from_which_chrom(i, true);
	const auto	probs_pat = this->probs_from_which_chrom(i, false);
	return this->likelihood_each(this->gt(i), probs_mat, probs_pat,
											mat_phasing, pat_phasing);
}

// 親のGTがひっくり返っているかひっくり返っていないのか仮定して尤度を計算する
// GTが0/1のとき、phasingが0なら0|1、1なら1|0
// GTが0|0のとき、どちらでも0|0
double VCFFillable::RecordSet::compute_phasing_likelihood(int mat_phasing,
														int pat_phasing) const {
	double	ll = 0.0;
	for(int i = 2; i < (int)record->num_samples(); ++i) {
		if(record->get_GT(i) == "./.")
			ll += log(0.0001);	// 本来要らないがPython版と合わせるため
		else
			ll += compute_phasing_likelihood_each(i, mat_phasing, pat_phasing);
	}
	return ll;
}

pair<int, int> VCFFillable::RecordSet::select_phasing(
									const vector<pair<int, int>>& candidates) {
	if(candidates.size() == 1)
		return candidates.front();
	
	const int	mat_int_gt = record->mat_int_gt();
	const int	pat_int_gt = record->pat_int_gt();
	
	vector<tuple<int, int, int>>	v;	// [(score, mat_phasing, pat_phasing)]
	for(auto p = candidates.begin(); p != candidates.end(); ++p) {
		const int	mat_int_gt1 = (p->first >> 1) + (p->first & 1);
		const int	pat_int_gt1 = (p->second >> 1) + (p->second & 1);
		const int	score = abs(mat_int_gt1 - mat_int_gt) +
							abs(pat_int_gt1 - pat_int_gt);
		v.push_back(make_tuple(score, p->first, p->second));
	}
	// これだと、同じスコアならphasingが小さい方から選んでいる
	// 本当はランダム的に選びたい
	const auto	p = std::min_element(v.begin(), v.end());
	return make_pair(get<1>(*p), get<2>(*p));
}

pair<int, int> VCFFillable::RecordSet::determine_phasing_core(
								const vector<tuple<double, int, int>>& lls) {
	// 最大に近いllを集める
	vector<pair<int, int>>	candidates;
	const double	max_ll = get<0>(lls.back());
	for(size_t i = lls.size() - 1; i < lls.size(); --i) {
		const double	ll = get<0>(lls[i]);
		if(ll > max_ll - 1e-9)
			candidates.push_back(make_pair(get<1>(lls[i]), get<2>(lls[i])));
		else
			break;
	}
	
	return this->select_phasing(candidates);
}

void VCFFillable::RecordSet::determine_phasing() {
	if(this->record == NULL)
		return;
	
	const vector<pair<int, int>>	phasing = record->possible_phasings();
//	double	max_ll = -std::numeric_limits<double>::max();
	vector<tuple<double, int, int>>	lls;
	for(auto p = phasing.begin(); p != phasing.end(); ++p) {
		const double	ll = compute_phasing_likelihood(p->first, p->second);
		lls.push_back(make_tuple(ll, p->first, p->second));
	}
	
	// 最大のllとほとんど同じllがあったとき、
	std::sort(lls.begin(), lls.end());
	const auto	p = this->determine_phasing_core(lls);
	const int	mat_phasing = p.first;
	const int	pat_phasing = p.second;
	
	// この処理は本当に合っているのか
	static const string	gts[] = { "0|0", "1|0", "0|1", "1|1" };
	record->set_mat_GT(gts[mat_phasing]);
	record->set_pat_GT(gts[pat_phasing]);
}

int VCFFillable::RecordSet::select_from(int from1, int from2,
										const VCFRecord *record1,
										const VCFRecord *record2) const {
	// どちらかは0でない前提
	if(from1 == 0)
		return from2;
	else if(from2 == 0)
		return from1;
	else {	// 両側にRecordがある
		// 近い方を選ぶ
		if(record->pos() * 2 < record1->pos() + record2->pos())
			return from1;
		else
			return from2;
	}
}

string VCFFillable::RecordSet::modify_gt(size_t i) {
	const int	prev_mat_from = this->prev_mat_from(i);
	const int	next_mat_from = this->next_mat_from(i);
	const int	prev_pat_from = this->prev_pat_from(i);
	const int	next_pat_from = this->next_pat_from(i);
	if((prev_mat_from == 0 && next_mat_from == 0) ||
							(prev_pat_from == 0 && next_pat_from == 0)) {
		return record->get_gt(i);
	}
	
	vector<pair<int,int>>	pairs_ = {
				pair<int,int>(prev_mat_from, prev_pat_from),
				pair<int,int>(prev_mat_from, next_pat_from),
				pair<int,int>(next_mat_from, prev_pat_from),
				pair<int,int>(next_mat_from, next_pat_from)
	};
	vector<pair<int,int>>	pairs__;
	for(auto p = pairs_.begin(); p != pairs_.end(); ++p) {
		if(p->first != 0 && p->second != 0)
			pairs__.push_back(*p);
	}
	const auto	pairs = Common::unique_vector(pairs__);
	const auto	p = select_pair(pairs, i);
	const int	mat_from = p.first;
	const int	pat_from = p.second;
	if(mat_from == 0 && pat_from == 0) {
		if(prev_mat_from != 0 || next_mat_from != 0) {
			const int	mat_from_selected = select_from(
											prev_mat_from, next_mat_from,
											prev_mat_record, next_mat_record);
			return record->gt_from_mat(mat_from_selected, i + 9);
		}
		else if(prev_pat_from != 0 || next_pat_from != 0) {
			const int	pat_from_selected = select_from(
											prev_pat_from, next_pat_from,
											prev_pat_record, next_pat_record);
			return record->gt_from_pat(pat_from_selected, i + 9);
		}
		else {
			// if both are 0, do nothing
			return record->get_GT(i);
		}
	}
	else {
		return record->gt_from_parent(mat_from, pat_from);
	}
}

void VCFFillable::RecordSet::impute_core() {
	STRVEC	new_gts;
	for(size_t i = 2; i < record->num_samples(); ++i) {
		new_gts.push_back(modify_gt(i));
	}
	record->modify_gts(new_gts);
	record->modify_parents_type();
}

int VCFFillable::RecordSet::select_mat(const vector<Pair>& pairs) const {
	if(pairs.size() == 1)
		return pairs.front().first;
	else if(this->prev_mat_record == NULL || this->next_mat_record == NULL)
		return 0;
	else if(this->is_mat_prev_near())
		return pairs.front().first;		// prev_mat_from
	else
		return pairs.back().first;		// next_mat_from
}

void VCFFillable::RecordSet::impute_NA_mat_each(size_t i) const {
	const int	prev_mat_from = VCFFillableRecord::from_which_chrom_mat(
													this->prev_mat_record, i);
	const int	next_mat_from = VCFFillableRecord::from_which_chrom_mat(
													this->next_mat_record, i);
	
	// patはホモだから1でも2でもよい
	vector<Pair>	pairs;	// [(mat_from, pat_from)]
	if(this->prev_mat_record != NULL)
		pairs.push_back(Pair(prev_mat_from, 1));
	if(this->next_mat_record != NULL && next_mat_from != prev_mat_from)
		pairs.push_back(Pair(next_mat_from, 1));
	if(pairs.empty())
		return;
	
	const int	mat_from = this->select_mat(pairs);
	this->record->set_GT(i, this->record->gt_from_parent(mat_from, 1));
}

int VCFFillable::RecordSet::select_pat(const vector<Pair>& pairs) const {
	if(pairs.size() == 1)
		return pairs.front().second;
	else if(this->prev_pat_record == NULL || this->next_pat_record == NULL)
		return 0;
	else if(this->is_pat_prev_near())
		return pairs.front().second;	// prev_pat_from
	else
		return pairs.back().second;		// next_pat_from
}

void VCFFillable::RecordSet::impute_NA_pat_each(size_t i) const {
	const int	prev_pat_from = VCFFillableRecord::from_which_chrom_pat(
													this->prev_pat_record, i);
	const int	next_pat_from = VCFFillableRecord::from_which_chrom_pat(
													this->next_pat_record, i);
	
	vector<Pair>	pairs;
	if(this->prev_pat_record != NULL)
		pairs.push_back(Pair(1, prev_pat_from));
	if(this->next_pat_record != NULL && next_pat_from != prev_pat_from)
		pairs.push_back(Pair(1, next_pat_from));
	if(pairs.empty())
		return;
	
	const int	pat_from = this->select_pat(pairs);
	this->record->set_GT(i, this->record->gt_from_parent(1, pat_from));
}


//////////////////// VCFFillable ////////////////////

VCFFillable::VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
								std::vector<VCFFillableRecord *> rs) :
				VCFBase(h, s), VCFSmallBase(), VCFFamilyBase(), records(rs) { }

VCFFillable::~VCFFillable() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFFillable::modify() {
	vector<Group>	groups = group_records();
	for(size_t i = 0U; i < groups.size(); ++i) {
		if(groups[i].second.front()->is_fillable_type())
			this->phase(i, groups, true);
	}
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
	
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->fill_PGT();
}

void VCFFillable::phase_hetero_hetero() {
	// typeが'IMPUTABLE', 'MAT', 'PAT', 'FIXED'でrecordを分ける
	vector<Group>	groups = this->group_records();
	for(size_t i = 0; i < groups.size(); ++i) {
		if(groups[i].first == FillType::IMPUTABLE)
			this->phase(i, groups, false);
	}
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
}

VCFFillable *VCFFillable::create_from_header() const {
	vector<VCFFillableRecord *>	rs;
	VCFFillable	*vcf = new VCFFillable(header, samples, rs);
	copy_chrs(vcf);
	return vcf;
}

vector<VCFFillable::Group> VCFFillable::group_records() const {
	vector<Group>	groups;
	auto	current_type = records.front()->get_type();
	vector<VCFFillableRecord *>	group(1U, records.front());
	for(auto p = records.begin() + 1; p != records.end(); ++p) {
		auto	*record = *p;
		if(record->get_type() != current_type) {
			groups.push_back(Group(current_type, group));
			group.clear();
			current_type = record->get_type();
		}
		group.push_back(record);
	}
	groups.push_back(Group(current_type, group));
	return groups;
}

VCFFillableRecord *VCFFillable::find_prev_record(FillType type, int i,
											const vector<Group>& groups) const {
	const string&	chr = groups[i].second.front()->chrom();
	for(int j = i - 1; j >= 0; --j) {
		if(groups[j].second.front()->chrom() != chr)
			return NULL;
		else if(groups[j].first == type)
			return groups[j].second.back();
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_next_record(FillType type, int i,
											const vector<Group>& groups) const {
	const string&	chr = groups[i].second.front()->chrom();
	for(int j = i + 1; j < (int)groups.size(); ++j) {
		if(groups[j].second.front()->chrom() != chr)
			return NULL;
		else if(groups[j].first == type)
			return groups[j].second.front();
	}
	return NULL;
}

void VCFFillable::phase(int i, const vector<Group>& groups,
										bool necessary_parents_phasing) {
	const auto	mat = FillType::MAT;
	const auto	pat = FillType::PAT;
	auto	*prev_mat_record = find_prev_record(mat, i, groups);
	auto	*next_mat_record = find_next_record(mat, i, groups);
	auto	*prev_pat_record = find_prev_record(pat, i, groups);
	auto	*next_pat_record = find_next_record(pat, i, groups);
	auto	group_records = groups[i].second;
	for(auto p = group_records.begin(); p != group_records.end(); ++p) {
		auto	*record = *p;
		RecordSet	record_set(record, prev_mat_record, next_mat_record,
											prev_pat_record, next_pat_record);
		if(necessary_parents_phasing)
			record_set.determine_phasing();
		record_set.impute_core();
	}
}

template<typename Iter>
VCFFillableRecord *VCFFillable::find_neighbor_same_type_record(
							size_t i, size_t c, Iter first, Iter last) const {
	const FillType	type = records[i]->get_type();
	const string&	chromosome = records[i]->chrom();
	for(auto p = first; p != last; ++p) {
		auto	*record = *p;
		if(record->chrom() != chromosome)
			return NULL;
		else if(record->get_type() == type && record->get_GT(c-9) != "./.")
			return record;
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_prev_same_type_record(
													size_t i, size_t c) const {
	if(i == 0U)
		return NULL;
	
	return find_neighbor_same_type_record(i, c, records.rend() - i + 1,
																records.rend());
}

VCFFillableRecord *VCFFillable::find_next_same_type_record(
													size_t i, size_t c) const {
	if(i == records.size() - 1)
		return NULL;
	
	const auto	first = records.begin() + i + 1;
	const auto	last = records.end();
	return find_neighbor_same_type_record(i, c, first, last);
}

const VCFFillable::RecordSet *VCFFillable::create_recordset(
										size_t i, size_t c, bool is_mat) const {
	auto	*record = records[i];
	auto	*prev_record = find_prev_same_type_record(i, c);
	auto	*next_record = find_next_same_type_record(i, c);
	if(is_mat)
		return new RecordSet(record, prev_record, next_record, NULL, NULL);
	else
		return new RecordSet(record, NULL, NULL, prev_record, next_record);
}

void VCFFillable::impute_NA_mat_each(size_t i, size_t c) {
	const RecordSet	*rs = this->create_recordset(i, c, true);
	rs->impute_NA_mat_each(c - 9);
	delete rs;
}

void VCFFillable::impute_NA_mat(size_t i) {
	auto	*record = this->records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_mat_each(i, c);
	}
}

void VCFFillable::impute_NA_pat_each(size_t i, size_t c) {
	const RecordSet	*rs = this->create_recordset(i, c, false);
	rs->impute_NA_pat_each(c - 9);
	delete rs;
}

void VCFFillable::impute_NA_pat(size_t i) {
	auto	*record = this->records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_pat_each(i, c);
	}
}

pair<int, int> VCFFillable::find_prev_mat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->mat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_mat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->mat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_prev_pat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->pat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_pat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->pat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

int VCFFillable::select_from(const pair<int, int>& f1,
							 const pair<int, int>& f2, int i) const {
	const int	i1 = f1.first;
	const int	from1 = f1.second;
	const int	i2 = f2.first;
	const int	from2 = f2.second;
	if(from1 == 0 && from2 == 0) {
		// 前後がないとき乱数的に決める
		const auto	*r0 = this->records[i];
		return r0->pos() % 2 + 1;
	}
	if(from1 == from2)
		return from1;
	else if(from2 == 0)
		return from1;
	else if(from1 == 0)
		return from2;
	else {
		// 最後は物理距離で決める
		const auto	*r0 = this->records[i];
		const auto	*r1 = this->records[i1];
		const auto	*r2 = this->records[i2];
		if(r0->pos() * 2 <= r1->pos() + r2->pos())
			return from1;
		else
			return from2;
	}
}

int VCFFillable::find_mat_from(int i, int c) const {
	return this->select_from(this->find_prev_mat_from(i, c),
							 this->find_next_mat_from(i, c), i);
}

int VCFFillable::find_pat_from(int i, int c) const {
	return this->select_from(this->find_prev_pat_from(i, c),
							 this->find_next_pat_from(i, c), i);
}

void VCFFillable::impute_others(int i) {
	auto	*record = this->records[i];
	const bool	mat_homo = record->is_homo(0);
	const bool	pat_homo = record->is_homo(1);
	for(size_t c = 11; c != record->get_v().size(); ++c) {
		if(record->get_v()[c].c_str()[1] != '/' &&
								record->get_int_gt(c-9) != -1)
			continue;
		const int	mat_from = mat_homo ? 1 : this->find_mat_from(i, c);
		const int	pat_from = pat_homo ? 1 : this->find_pat_from(i, c);
		record->set_GT(c-9, record->gt_from_parent(mat_from, pat_from));
	}
}

VCFFillable *VCFFillable::fill(const vector<VCFHeteroHomo *>& vcfs,
					const vector<VCFImpFamilyRecord *>& records, bool all_out) {
	vector<VCFFillableRecord *>	merged_records;
	if(all_out)
		merged_records = VCFFillable::merge_records(vcfs, records, all_out);
	else {
		vector<VCFImpFamilyRecord *>	records_;
		for(auto p = records.begin(); p != records.end(); ++p) {
			if(!(*p)->is_fixed())
				records_.push_back(*p);
		}
		merged_records = merge_records(vcfs, records_, all_out);
	}
	VCFFillable	*vcf = new VCFFillable(vcfs.front()->get_header(),
								vcfs.front()->get_samples(), merged_records);
	vcf->modify();
	return vcf;
}

vector<VCFFillableRecord *> VCFFillable::merge_records(
									const vector<VCFHeteroHomo *>& vcfs,
									const vector<VCFImpFamilyRecord *>& records,
									bool all_out) {
	vector<VCFFillableRecord *>	all_records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		auto	*vcf = *p;
		auto&	hh_records = vcf->get_records();
		for(auto q = hh_records.begin(); q != hh_records.end(); ++q)
			all_records.push_back(VCFFillableRecord::convert(*q));
	}
	
	for(auto p = records.begin(); p != records.end(); ++p) {
		if(all_out || (*p)->is_fillable())
			all_records.push_back(VCFFillableRecord::convert(*p));
	}
	std::sort(all_records.begin(), all_records.end(),
				[](const VCFFillableRecord *r1, const VCFFillableRecord *r2) {
					return r1->get_index() < r2->get_index(); });
	return all_records;
}

vector<vector<VCFFillableRecord *>> VCFFillable::collect_records(
										const vector<VCFFillable *>& vcfs,
										bool all_out) {
	vector<vector<VCFFillableRecord *>>	rss;
	if(all_out) {
		for(size_t i = 0; i < vcfs.front()->size(); ++i) {
			vector<VCFFillableRecord *>	rs;
			for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
				rs.push_back((*p)->get_fillable_record(i));
			rss.push_back(rs);
		}
	}
	else {
		int	max_index = 0;
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
			max_index = max(max_index, (*p)->get_records().back()->get_index());
		}
		
		vector<size_t>	js(vcfs.size());	// 各VCFの現在の位置
		for(int i = 0; i <= max_index; ++i) {
			// このindexで全ての家系でrecordが揃っていればrssに登録する
			vector<VCFFillableRecord *>	rs(vcfs.size(), NULL);
			for(size_t k = 0; k < vcfs.size(); ++k) {
				VCFFillable	*vcf = vcfs[k];
				if(js[k] == vcf->size())
					continue;
				VCFFillableRecord	*record = vcf->get_fillable_record(js[k]);
				if(record->get_index() == i) {
					rs[k] = record;
					js[k] += 1;
				}
				
				if(std::all_of(rs.begin(), rs.end(),
									[](auto *r) { return r != NULL; }))
					rss.push_back(rs);
			}
		}
	}
	return rss;
}

pair<STRVEC, vector<vector<pair<int, int>>>>
VCFFillable::integrate_samples(const vector<STRVEC>& sss,
									const STRVEC& orig_samples) {
	// { sample: (index of family, index of inner family) }
	map<string, vector<pair<int, int>>>	dic;
	for(size_t i = 0; i < sss.size(); ++i) {
		const vector<string>&	samples = sss[i];
		for(size_t j = 0; j < samples.size(); ++j)
			dic[samples[j]].push_back(P(i, j));
	}
	
	STRVEC	new_samples;
	vector<vector<pair<int, int>>>	pos_samples;
	for(auto p = orig_samples.begin(); p != orig_samples.end(); ++p) {
		const string&	sample = *p;
		if(dic.find(sample) != dic.end()) {
			new_samples.push_back(sample);
			pos_samples.push_back(dic[sample]);
		}
	}
	
	return P(new_samples, pos_samples);
}

// 重複したサンプルが一つになるようにVCFを統合する
VCFSmall *VCFFillable::integrate(const VCFFillable *vcf,
								const vector<vector<VCFFillableRecord *>>& rss,
								const STRVEC& orig_samples) {
	vector<STRVEC>	sss;
	for(auto p = rss.front().begin(); p != rss.front().end(); ++p) {
		const STRVEC&	samples = (*p)->get_samples();
		sss.push_back(samples);
	}
	
	const auto	q = integrate_samples(sss, orig_samples);
	const STRVEC&	new_samples = q.first;
	// サンプルが各Familyの何番目にあるか
	const auto&		pos_samples = q.second;
	
	const vector<STRVEC>	header = vcf->trim_header(new_samples);
	VCFSmall	*new_vcf = new VCFSmall(header, new_samples,
											vector<VCFRecord *>());
	const STRVEC&	samples = new_vcf->get_samples();
	vector<VCFRecord *>	records;
	for(auto p = rss.begin(); p != rss.end(); ++p) {
		auto	*r = VCFFillableRecord::integrate(*p, samples, pos_samples);
		records.push_back(r);
	}
	new_vcf->add_records(records);
	return new_vcf;
}

VCFSmall *VCFFillable::merge(const vector<VCFFillable *>& vcfs,
											const STRVEC& orig_samples) {
	const auto	rss = collect_records(vcfs, true);
	return integrate(vcfs.front(), rss, orig_samples);
}

void VCFFillable::fill_in_thread(void *config) {
	const auto	*c = (ConfigFillThread *)config;
	const size_t	n = c->size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	vcfs = c->items[i].first;
		auto	records = c->items[i].second;
		auto	result = VCFFillable::fill(vcfs, records, c->all_out);
		c->filled_vcfs[i] = result;
	}
}

vector<VCFFillable *> VCFFillable::fill_parellel(vector<Item>& items,
														int num_threads) {
	vector<VCFFillable *>	results(items.size());
	
	const int	T = min((int)items.size(), num_threads);
	vector<ConfigFillThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigFillThread(items, true, i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&fill_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		fill_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return results;
}

void VCFFillable::delete_items(const vector<Item>& items) {
	for(auto q = items.begin(); q != items.end(); ++q) {
		for(auto r = q->first.begin(); r != q->first.end(); ++r)
			delete *r;
		for(auto r = q->second.begin(); r != q->second.end(); ++r)
			delete *r;
	}
}

vector<VCFFillable *> VCFFillable::fill_all(
							map<Parents, vector<VCFHeteroHomo *>>& imputed_vcfs,
							ImpRecords& other_records, int num_threads) {
	vector<Item>	items;
	for(auto q = imputed_vcfs.begin(); q != imputed_vcfs.end(); ++q) {
		const Parents&	parents = q->first;
		vector<VCFHeteroHomo *>&	vcfs = q->second;
		items.push_back(make_pair(vcfs, other_records[parents]));
	}
	vector<VCFFillable *>	filled_vcfs = fill_parellel(items, num_threads);
	delete_items(items);
	return filled_vcfs;
}
