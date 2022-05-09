#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "VCFFillable.h"
#include "common.h"

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

bool Genotype::is_consistent(const Genotype& mat_gt, const Genotype& pat_gt,
												bool considers_phasing) const {
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

bool Genotype::is_valid(const string& gt) {
	return (gt.c_str()[0] == '0' || gt.c_str()[0] == '1') &&
			(gt.c_str()[2] == '0' || gt.c_str()[2] == '1');
}

int Genotype::sum_gt(const string& gt) {
	return (int)((gt.c_str()[0] - '0') + (gt.c_str()[2] - '0'));
}


//////////////////// VCFFillableRecord ////////////////////

STRVEC VCFFillableRecord::prog_gts() const {
	STRVEC	vec(v.size() - 11);
	std::copy(v.begin() + 11, v.end(), vec.begin());
	return vec;
}

VCFFillableRecord *VCFFillableRecord::copy() const {
	return new VCFFillableRecord(v, samples, type);
}

string VCFFillableRecord::gt_from_parent(int mat_from, int pat_from) const {
	stringstream	ss;
	ss << v[9].c_str()[mat_from*2-2] << '|' << v[10].c_str()[pat_from*2-2];
	return ss.str();
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

int VCFFillableRecord::segregation_type() const {
	const auto	pss = make_probability_table();
	const auto	gts = progeny_gts();	// [gt: int(0|1|2)]
	vector<int>	ns(3, 0);
	for(auto p = gts.begin(); p != gts.end(); ++p) {
		if(0 <= *p && *p < 3)
			ns[*p] += 1;
	}
	
	int	seg_type = 0;	// temporary
	double	max_lls = -1e+308;
	for(int i = 0; i < 6; ++i) {
		double	lls = 0.0;
		for(int j = 0; j < 3; ++j)
			lls += ns[j] * log(pss[i][j]);
		if(lls > max_lls) {
			seg_type = i;
			max_lls = lls;
		}
	}
	return seg_type;
}

vector<size_t> VCFFillableRecord::conflicted_progeny_indices() const {
	const auto	vec_gts = this->gts();
	const Genotype	mat_gt(vec_gts[0]);
	const Genotype	pat_gt(vec_gts[1]);
	
	vector<size_t>	indices;
	for(auto p = vec_gts.begin() + 2; p != vec_gts.end(); ++p) {
		const Genotype	prog_gt(*p);
		// phasingを考慮に入れずに矛盾をカウント
		// phasingはあとで修正する
		if(!prog_gt.is_consistent(mat_gt, pat_gt, false))
			indices.push_back(p - vec_gts.begin() - 2);
	}
	return indices;
}

bool VCFFillableRecord::modify_parents() {
	map<string,int>	counter;
	const auto	vec_gts = this->gts();
	for(auto p = vec_gts.begin(); p != vec_gts.end(); ++p)
		counter[p->substr(0, 3)] += 1;
	
	const int	c00 = counter["0/0"];
	const int	c01 = counter["0/1"];
	const int	c11 = counter["1/1"];
	if(c00 * 10 >= (c00 + c01 + c11) * 9) {
		this->set_mat_GT("0/0");
		this->set_pat_GT("0/0");
		return true;
	}
	else if(c01 * 10 >= (c00 + c01 + c11) * 9) {
		// 親は0/0と1/1のはず
		if(this->get_GT(9) == "0/0") {
			if(this->get_GT(10) == "0/0") {
				return false;	// どちらを修正していいのか分からない
			}
			else {
				this->set_pat_GT("1/1");
				return true;
			}
		}
		else if(this->get_GT(9) == "1/1") {
			if(this->get_GT(10) == "1/1") {
				return false;
			}
			else {
				this->set_pat_GT("0/0");
				return true;
			}
		}
		else {
			if(this->get_GT(10) == "0/0") {
				this->set_mat_GT("1/1");
				return true;
			}
			else if(this->get_GT(10) == "1/1") {
				this->set_mat_GT("0/0");
				return true;
			}
			else {
				return false;
			}
		}
	}
	else if(c11 * 10 >= (c00 + c01 + c11) * 9) {
		this->set_mat_GT("1/1");
		this->set_pat_GT("1/1");
		return true;
	}
	else {
		return false;
	}
}

bool VCFFillableRecord::is_valid() {
	const int	gt_mat = this->get_int_gt(0);
	const int	gt_pat = this->get_int_gt(1);
	
	// 親に対して矛盾した子が10%以上
	if(conflicted_progeny_indices().size()*10 >= num_progenies())
		return modify_parents();
	
	// 親のどちらか一方がヘテロ
	if(this->is_homo(0) ^ this->is_homo(1))
		return false;
	
	// 尤度計算して両親と合っているか
	const int	s = segregation_type();
	const int	g1 = min(gt_mat, gt_pat);
	const int	g2 = max(gt_mat, gt_pat);
	if(g1*(5-g1)/2 + g2 != s)
		return false;
	
	return true;
}

void VCFFillableRecord::modify() {
	if(this->type != RecordType::FILLED)
		return;
	else if(!this->is_valid())
		// あとで各家系を統合するときに楽なので、
		// 消すのではなく無効なレコードにする
		this->disable();
	else if(this->is_homo(0) && this->is_homo(1))
		this->phase();
}

void VCFFillableRecord::phase() {
	// 両親ともホモになっている前提
	const string	gt_mat = this->v[9].substr(0, 1);
	const string	gt_pat = this->v[10].substr(0, 1);
	this->set_mat_GT(gt_mat + "|" + gt_mat);
	this->set_pat_GT(gt_pat + "|" + gt_pat);
	for(size_t i = 11U; i < this->v.size(); ++i)
		this->set_GT(i, gt_mat + "|" + gt_pat);
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
		set_GT(9, v[9].substr(2, 1) + "|" + v[9].substr(0, 1));
	if(inv_pat)
		set_GT(10, v[10].substr(2, 1) + "|" + v[10].substr(0, 1));
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
		if(gts[i] != "0/1")
			num += 1;
		if(!is_same_gts(v[i+11], gts[i]))
			dist += 1;
	}
	return dist < num / 2;
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

int VCFFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	int j = is_mat ? 0 : 1;
	const string&	parent_gt = this->get_gt(j);
	return parent_gt.c_str()[0] == this->get_gt(i).c_str()[j*2] ? 1 : 2;
}

VCFFillableRecord *VCFFillableRecord::from_VCFFamilyRecord(
											const VCFFamilyRecord *record) {
	auto	type = record->is_homo(1) ? RecordType::MAT : RecordType::PAT;
	return new VCFFillableRecord(record->get_v(), record->get_samples(), type);
}

void VCFFillableRecord::set(const STRVEC& new_v, RecordType new_type) {
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

VCFFillable::Pair VCFFillable::RecordSet::select_nearest_froms(
								const vector<Pair>& pairs, size_t i) const {
	if(pairs.size() == 4U) {
		return Pair(is_mat_prev_near() ? prev_mat_from(i) : next_mat_from(i),
					is_pat_prev_near() ? prev_pat_from(i) : next_pat_from(i));
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
	else {	// ねじれている
		if(next_pat_record->pos() - prev_mat_record->pos() <
								next_mat_record->pos() - prev_pat_record->pos())
			return Pair(prev_mat_from(i), next_pat_from(i));
		else
			return Pair(next_mat_from(i), prev_pat_from(i));
	}
}

VCFFillable::Pair VCFFillable::RecordSet::select_pair(const vector<Pair>& pairs,
												size_t i, bool selected) const {
	if(pairs.empty())
		return Pair(0, 0);
	else if(pairs.size() == 1U)
		return pairs.front();
	else if(Genotype::is_valid(record->get_gt(i)))
		return select_nearest_froms(pairs, i);
	else if(selected)
		return select_nearest_froms(pairs, i);
	
	vector<Pair>	new_pairs;
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


//////////////////// VCFFillable ////////////////////

vector<VCFRecord *> VCFFillable::to_VCFRecord(vector<VCFFillableRecord *>& rs) {
	vector<VCFRecord *>	records(rs.size());
	std::copy(rs.begin(), rs.end(), records.begin());
	return records;
}

VCFFillable *VCFFillable::insert_positions(const vector<Position>& positions) {
	// とりあえず、ポジションだけの空のレコードを作る
	vector<VCFFillableRecord *>	new_records;
	size_t	i = 0U;
	auto	pos1 = record_position(*records[i]);
	for(auto p = positions.begin(); p != positions.end(); ++p) {
		const auto&	pos2 = *p;
		if(pos1.first == get<0>(pos2) && pos1.second == get<1>(pos2)) {
			new_records.push_back(records[i]->copy());
			// proceed
			i += 1;
			if(i < records.size())
				pos1 = record_position(*records[i]);
		}
		else {
			STRVEC	vec(9);
			vec[0] = this->chr(get<0>(pos2));
			stringstream	ss2;
			ss2 << get<1>(pos2);
			vec[1] = ss2.str();
			new_records.push_back(new VCFFillableRecord(vec, samples,
										VCFFillableRecord::RecordType::FILLED));
		}
	}
	
	return new VCFFillable(header, samples, new_records);
}

vector<VCFFillable::Group> VCFFillable::group_records() const {
	vector<Group>	groups;
	VCFFillableRecord::RecordType	current_type = records.front()->get_type();
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

VCFFillableRecord *VCFFillable::find_prev_record(
								VCFFillableRecord::RecordType type,
								int i, const vector<Group>& groups) const {
	for(int j = i - 1; j >= 0; --j) {
		if(groups[j].first == type)
			return groups[j].second.back();
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_next_record(
								VCFFillableRecord::RecordType type,
								int i, const vector<Group>& groups) const {
	for(int j = i + 1; j < (int)groups.size(); ++j) {
		if(groups[j].first == type)
			return groups[j].second.front();
	}
	return NULL;
}

// phasingされている前提
int VCFFillable::from_which_chrom(const VCFFillableRecord *record,
											size_t i, bool is_mat) const {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}

vector<double> VCFFillable::probs_from_which_chrom(
									int prev_chrom, int next_chrom) const {
	double	prob1;
	if(prev_chrom == 1 && next_chrom == 1)
		prob1 = 0.99;
	else if(prev_chrom == 2 && next_chrom == 2)
		prob1 = 0.01;
	else
		prob1 = 0.5;
	
	vector<double>	probs = { prob1, 1.0 - prob1 };
	return probs;
}

double VCFFillable::likelihood_each(const string& gt,
									const vector<double>& probs_mat,
									const vector<double>& probs_pat,
									int mat_phasing, int pat_phasing) const {
	const int	sum = Genotype::sum_gt(gt);
	double	likelihood = 0.0;
	for(int k = 0; k < 4; ++k) {
		const int	i = k >> 1;
		const int	j = k & 1;
		if((((i + mat_phasing) & 1) + ((j + pat_phasing) & 1)) == sum)
			likelihood += probs_mat[i] * probs_pat[j];
	}
	return log(likelihood);
}

// mat_phasing : '0|1' or '1|0'
double VCFFillable::compute_phasing_likelihood_each(RecordSet& rs, size_t i,
									int mat_phasing, int pat_phasing) const {
	const auto	probs_mat = probs_from_which_chrom(rs.prev_mat_from(i),
													rs.next_mat_from(i));
	const auto	probs_pat = probs_from_which_chrom(rs.prev_pat_from(i),
													rs.next_pat_from(i));
	return likelihood_each(rs.gt(i), probs_mat, probs_pat,
									mat_phasing, pat_phasing);
}

// mat_phasing : '0|1' or '1|0'
double VCFFillable::compute_phasing_likelihood(RecordSet& rs, int mat_phasing,
														int pat_phasing) const {
	double	ll = 0.0;
	for(int i = 2; i < (int)samples.size(); ++i) {
		if(rs.record->get_GT(i) == "./.")
			continue;
		ll += compute_phasing_likelihood_each(rs, i, mat_phasing, pat_phasing);
	}
	return ll;
}

void VCFFillable::determine_phasing(RecordSet& rs) {
	int	mat_phasing = 0;
	int	pat_phasing = 0;
	double	max_ll = compute_phasing_likelihood(rs, 0, 0);
	for(int i = 1; i < 4; ++i) {
		const double	ll = compute_phasing_likelihood(rs, i >> 1, i & 1);
		if(ll > max_ll) {
			mat_phasing = i >> 1;
			pat_phasing = i & 1;
			max_ll = ll;
		}
	}
	
	rs.record->set_mat_GT(mat_phasing == 0 ? "0|1" : "1|0");
	rs.record->set_pat_GT(pat_phasing == 0 ? "0|1" : "1|0");
}

string VCFFillable::modify_gt(RecordSet& rs, size_t i) {
	const int	prev_mat_from = rs.prev_mat_from(i);
	const int	next_mat_from = rs.next_mat_from(i);
	const int	prev_pat_from = rs.prev_pat_from(i);
	const int	next_pat_from = rs.next_pat_from(i);
	if((prev_mat_from == 0 && next_mat_from == 0) ||
							(prev_pat_from == 0 && next_pat_from == 0)) {
		return rs.record->get_gt(i);
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
	const auto	p = rs.select_pair(pairs, i);
	return rs.record->gt_from_parent(p.first, p.second);
}

void VCFFillable::impute_core(RecordSet& rs) {
	STRVEC	new_gts;
	for(size_t i = 2U; i < samples.size(); ++i) {
		new_gts.push_back(modify_gt(rs, i));
	}
	rs.record->modify_gts(new_gts);
}

void VCFFillable::phase(int i, const vector<Group>& groups) {
	const auto	mat = VCFFillableRecord::RecordType::MAT;
	const auto	pat = VCFFillableRecord::RecordType::PAT;
	auto	*prev_mat_record = find_prev_record(mat, i, groups);
	auto	*next_mat_record = find_next_record(mat, i, groups);
	auto	*prev_pat_record = find_prev_record(pat, i, groups);
	auto	*next_pat_record = find_next_record(pat, i, groups);
	auto	group_records = groups[i].second;
	for(auto p = group_records.begin(); p != group_records.end(); ++p) {
		auto	*record = *p;
		if(!record->is_homo(0) && !record->is_homo(1)) {
			RecordSet	record_set(record, prev_mat_record, next_mat_record,
											prev_pat_record, next_pat_record);
			determine_phasing(record_set);
			impute_core(record_set);
		}
	}
}

VCFFillableRecord *VCFFillable::find_prev_same_type_record(size_t i,
															size_t c) const {
	if(i == 0U)
		return NULL;
	
	const VCFFillableRecord::RecordType	type = records[i]->get_type();
	const string&	chromosome = records[i]->chrom();
	for(auto p = records.rend() - i + 1; p != records.rend(); ++p) {
		auto	*record = *p;
		if(record->chrom() != chromosome)
			return NULL;
		else if(record->get_type() == type && record->get_GT(c-9) != "./.")
			return record;
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_next_same_type_record(size_t i,
															size_t c) const {
	if(i == records.size() - 1)
		return NULL;
	
	const VCFFillableRecord::RecordType	type = records[i]->get_type();
	const string&	chromosome = records[i]->chrom();
	for(auto p = records.begin() + i + 1; p != records.end(); ++p) {
		auto	*record = *p;
		if(record->chrom() != chromosome)
			return NULL;
		else if(record->get_type() == type && record->get_GT(c-9) != "./.")
			return record;
	}
	return NULL;
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

int VCFFillable::select_mat(const vector<Pair>& pairs,
										const RecordSet *rs) const {
	if(pairs.size() == 1)
		return pairs.front().first;
	else if(rs->is_mat_prev_near())
		return pairs.front().first;		// prev_mat_from
	else
		return pairs.back().first;		// next_mat_from
}

void VCFFillable::impute_NA_mat_each(size_t i, size_t c) {
	const RecordSet	*rs = create_recordset(i, c, true);
	const int	prev_mat_from = from_which_chrom_mat(rs->prev_mat_record, c-9);
	const int	next_mat_from = from_which_chrom_mat(rs->next_mat_record, c-9);
	
	// なぜpatは1なのか
	vector<Pair>	pairs;
	if(prev_mat_from != 0)
		pairs.push_back(Pair(prev_mat_from, 1));
	if(next_mat_from != 0 && next_mat_from != prev_mat_from)
		pairs.push_back(Pair(next_mat_from, 1));
	if(pairs.empty())
		return;
	
	const int	mat_from = select_mat(pairs, rs);
	rs->record->set_GT(c, rs->record->gt_from_parent(mat_from, 1));
	delete rs;
}

void VCFFillable::impute_NA_mat(size_t i) {
	auto	*record = records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_mat_each(i, c);
	}
}

int VCFFillable::select_pat(const vector<Pair>& pairs,
										const RecordSet *rs) const {
	if(pairs.size() == 1)
		return pairs.front().second;
	else if(rs->is_pat_prev_near())
		return pairs.front().second;	// prev_pat_from
	else
		return pairs.back().second;		// next_pat_from
}

void VCFFillable::impute_NA_pat_each(size_t i, size_t c) {
	const RecordSet	*rs = create_recordset(i, c, false);
	const int	prev_pat_from = from_which_chrom_pat(rs->prev_pat_record, c-9);
	const int	next_pat_from = from_which_chrom_pat(rs->next_pat_record, c-9);
	
	// なぜpatは1なのか
	vector<Pair>	pairs;
	if(prev_pat_from != 0)
		pairs.push_back(Pair(1, prev_pat_from));
	if(next_pat_from != 0 && next_pat_from != prev_pat_from)
		pairs.push_back(Pair(1, next_pat_from));
	if(pairs.empty())
		return;
	
	const int	pat_from = select_pat(pairs, rs);
	rs->record->set_GT(c, rs->record->gt_from_parent(1, pat_from));
	delete rs;
}

void VCFFillable::impute_NA_pat(size_t i) {
	auto	*record = records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_pat_each(i, c);
	}
}

void VCFFillable::impute() {
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->modify();
	
	vector<Group>	groups = group_records();
	for(size_t i = 0U; i < groups.size(); ++i) {
		if(groups[i].second.front()->is_filled_type())
			this->phase(i, groups);
	}
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			impute_NA_mat(i);
		else if(record->is_pat_type())
			impute_NA_pat(i);
	}
	
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->fill_PGT();
}

VCFFillable *VCFFillable::convert(const VCFFamily *vcf) {
	vector<VCFFillableRecord *>	new_records;
	for(size_t i = 0U; i < vcf->size(); ++i) {
		const VCFFamilyRecord	*record = vcf->get_record(i);
		auto	new_record = VCFFillableRecord::from_VCFFamilyRecord(record);
		new_records.push_back(new_record);
	}
	return new VCFFillable(vcf->get_header(), vcf->get_samples(), new_records);
}

void VCFFillable::replace_filled_records(const vector<VCFFillable *>& vcfs,
															VCFHuge *orig_vcf) {
	POSITION	orig_pos(0, 0LL);
	VCFRecord	*orig_record = NULL;	// temporary
	VCFFillable	*vcf = vcfs.front();
	for(size_t i = 0U; i < vcf->size(); ++i) {
		const POSITION	pos = vcf->record_position(*(vcf->get_record(i)));
		if(pos > orig_pos) {
			// orig_posは最初は(0, 0)なので、orig_recordが不定なことはない
			orig_record = orig_vcf->proceed(pos);
			orig_pos = pos;
		}
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
			VCFFillableRecord	*record = (*p)->get_record(i);
			if(record->get_v().size() == 9U) {
				auto	v = orig_record->extract_v(record->get_samples());
				record->set(v, record->get_type());
			}
		}
	}
}

vector<VCFFillable::PosWithChr> VCFFillable::merge_positions(
								vector<VCFFillable *>::const_iterator first,
								vector<VCFFillable *>::const_iterator last) {
	vector<PosWithChr>	positions;
	if(last - first == 1) {
		auto	*vcf = *first;
		for(size_t i = 0U; i < vcf->size(); ++i) {
			auto	*record = vcf->get_record(i);
			auto	pos = vcf->record_position(*record);
			PosWithChr	position(pos.first, pos.second, record->chrom());
			positions.push_back(position);
		}
	}
	else {
		auto	mid = first + (last - first) / 2;
		const auto	positions1 = merge_positions(first, mid);
		const auto	positions2 = merge_positions(mid, last);
		size_t	k = 0U;
		size_t	l = 0U;
		while(k < positions1.size() && l < positions2.size()) {
			const PosWithChr	pos1 = positions1[k];
			const PosWithChr	pos2 = positions2[l];
			if(pos1 == pos2) {
				positions.push_back(pos1);
				k += 1;
				l += 1;
			}
			else if(pos1 < pos2) {
				positions.push_back(pos1);
				k += 1;
			}
			else {
				positions.push_back(pos2);
				l += 1;
			}
		}
		
		for( ; k < positions1.size(); ++k)
			positions.push_back(positions1[k]);
		for( ; l < positions2.size(); ++l)
			positions.push_back(positions2[l]);
	}
	return positions;
}

std::vector<VCFFillable::PosWithChr> VCFFillable::merge_positions(
									const std::vector<VCFFillable *>& vcfs) {
	return merge_positions(vcfs.begin(), vcfs.end());
}

STRVEC VCFFillable::join_samples(const vector<VCFFillable *>& vcfs) {
	// 本当は同じサンプルは排除したい
	STRVEC	samples;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const STRVEC&	vcf_samples = (*p)->get_samples();
		samples.insert(samples.end(), vcf_samples.begin(), vcf_samples.end());
	}
	return samples;
}

vector<STRVEC> VCFFillable::join_header(const vector<VCFFillable *>& vcfs,
													const STRVEC& samples) {
	const vector<STRVEC>&	vcf_header = vcfs.front()->get_header();
	vector<STRVEC>	header(vcf_header.begin(), vcf_header.end() - 1);
	STRVEC	bottom(vcf_header.back().begin(), vcf_header.back().begin() + 9);
	bottom.insert(bottom.end(), samples.begin(), samples.end());
	header.push_back(bottom);
	return header;
}

VCFSmall *VCFFillable::merge_vcfs(const vector<VCFFillable *>& vcfs) {
	const STRVEC	samples = join_samples(vcfs);
	const vector<STRVEC>	header = join_header(vcfs, samples);
	
	vector<VCFRecord *>	new_records;
	for(size_t i = 0U; i < vcfs.front()->size(); ++i) {
		vector<VCFFillableRecord *>	records;
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
			records.push_back((*p)->get_record(i));
		VCFRecord	*record = VCFFillableRecord::integrate_records(records);
		if(record != NULL)
			new_records.push_back(record);
	}
	return new VCFSmall(header, samples, new_records);
}

VCFSmall *VCFFillable::join_vcfs(const vector<VCFFamily *>& vcfs_,
											const string& path_VCF) {
	vector<VCFFillable *>	vcfs;
	for(auto p = vcfs_.begin(); p != vcfs_.end(); ++p)
		vcfs.push_back(VCFFillable::convert(*p));
	
	const vector<PosWithChr>	positions = merge_positions(vcfs);
	
	vector<VCFFillable *>	filled_vcfs;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFillable	*vcf = (*p)->insert_positions(positions);
		filled_vcfs.push_back(vcf);
	}
	Common::delete_all(vcfs);
	
	VCFHuge	*orig_vcf = VCFHuge::read(path_VCF);
	VCFFillable::replace_filled_records(filled_vcfs, orig_vcf);
	for(auto p = filled_vcfs.begin(); p != filled_vcfs.end(); ++p) {
		(*p)->impute();
	}
	
	return merge_vcfs(filled_vcfs);
}
