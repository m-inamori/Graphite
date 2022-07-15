#include <sstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <cassert>
#include "VCFFillable.h"
#include "TypeDeterminer.h"
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

tuple<int,int,int> VCFFillableRecord::count_gts() const {
	const auto	gts = get_int_gts();
	int	counter[3] = { 0, 0, 0 };
	for(auto p = gts.begin(); p != gts.end(); ++p)
		counter[*p] += 1;
	return make_tuple(counter[0], counter[1], counter[2]);
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

void VCFFillableRecord::phase() {
	// 両親ともホモになっている前提
	const string	gt_mat = this->v[9].substr(0, 1);
	const string	gt_pat = this->v[10].substr(0, 1);
	this->set_mat_GT(gt_mat + "|" + gt_mat);
	this->set_pat_GT(gt_pat + "|" + gt_pat);
	for(size_t i = 11U; i < this->v.size(); ++i)
		this->set_GT(i, gt_mat + "|" + gt_pat);
}

void VCFFillableRecord::modify() {
	if(this->type != RecordType::FILLED)
		return;
	else if(this->is_homo(0) && this->is_homo(1))
		this->phase();
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
		if(v[i+11] != "0/1")
			num += 1;
		if(!is_same_gts(gts[i], v[i+11]))
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
											const VCFFamilyRecord *record,
											const STRVEC& samples) {
	auto	type = record->is_homo(1) ? RecordType::MAT : RecordType::PAT;
	return new VCFFillableRecord(record->get_v(), samples, type);
}

void VCFFillableRecord::set(const STRVEC& new_v, RecordType new_type) {
	v = new_v;
	type = new_type;
	VCFFamilyRecord::set(new_v);
}

void VCFFillableRecord::determine_parents_type(int p) {
	if(p == -1) {
		disable();
		return;
	}
	
	const int	gt1 = p >> 2;
	const int	gt2 = p & 3;
	if(gt1 == gt2) {
		set_mat_int_GT(gt1);
		set_pat_int_GT(gt2);
		return;
	}
	
	// gt1とgt2を両親のどちらかに当てはめたときの当てはまり具合を調べる
	// 0: 当てはまらない 1: 当てはまる
	const int	mat_gt = mat_int_gt();
	const int	pat_gt = pat_int_gt();
	int	match1 = (mat_gt != gt1 ? 0 : 1) | (pat_gt != gt2 ? 0 : 2);
	int	match2 = (mat_gt != gt2 ? 0 : 1) | (pat_gt != gt1 ? 0 : 2);
	if(match1 == 3 || match2 == 3) {	// 合っている
		return;
	}
	else if(match1 != 0 && match2 == 0) {	// (gt1, gt2)はどちらかが合っている
		if(match1 == 1)		// matだけ合っている
			set_pat_int_GT(gt2);
		else				// patだけ合っている
			set_mat_int_GT(gt1);
	}
	else if(match1 == 0 && match2 != 0) {
		if(match2 == 1)		// matだけ合っている
			set_pat_int_GT(gt1);
		else
			set_mat_int_GT(gt2);
	}
	else {
		disable();
	}
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
	else if(!Genotype::is_valid(record->get_gt(i)))
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

VCFFillable::VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFillableRecord *> rs) :
					VCFSmall(h, s, vector<VCFRecord *>(rs.begin(), rs.end())),
					fillable_records(rs) { }

VCFFillable *VCFFillable::create_from_header() const {
	vector<VCFFillableRecord *>	rs;
	VCFFillable	*vcf = new VCFFillable(header, samples, rs);
	copy_chrs(vcf);
	return vcf;
}

void VCFFillable::set_records(const vector<VCFFillableRecord *>& rs) {
	fillable_records = rs;
	set_records_base(rs);
}

// 親に書くと子クラスのvectorを全て並べることになるので子に書く
void VCFFillable::set_records_base(const vector<VCFFillableRecord *>& rs) {
	records.insert(records.end(), rs.begin(), rs.end());
}

VCFFillable *VCFFillable::insert_positions(const vector<Position>& positions) {
	// 空のVCFFillableを作って、samplesを参照させる
	VCFFillable	*vcf = create_from_header();
	const STRVEC&	samples_ = vcf->get_samples();
	
	// とりあえず、ポジションだけの空のレコードを作る
	vector<VCFFillableRecord *>	new_records;
	size_t	i = 0U;
	auto	pos1 = record_position(*fillable_records[i]);
	for(auto p = positions.begin(); p != positions.end(); ++p) {
		const auto&	pos2 = *p;
		if(pos1.first == get<0>(pos2) && pos1.second == get<1>(pos2)
											&& i < fillable_records.size()) {
			// 新しいVCFを作る時はRecordも新しく作る
			new_records.push_back(fillable_records[i]->copy());
			// proceed
			i += 1;
			if(i < fillable_records.size())
				pos1 = record_position(*fillable_records[i]);
		}
		else {
			STRVEC	vec(9);
			vec[0] = this->chr(get<0>(pos2));
			stringstream	ss2;
			ss2 << get<1>(pos2);
			vec[1] = ss2.str();
			new_records.push_back(new VCFFillableRecord(vec, samples_,
										VCFFillableRecord::RecordType::FILLED));
		}
	}
	
	vcf->set_records(new_records);
	return vcf;
}

vector<VCFFillable::Group> VCFFillable::group_records() const {
	vector<Group>	groups;
	auto	current_type = fillable_records.front()->get_type();
	vector<VCFFillableRecord *>	group(1U, fillable_records.front());
	for(auto p = fillable_records.begin() + 1;
							p != fillable_records.end(); ++p) {
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
	const string&	chr = groups[i].second.front()->chrom();
	for(int j = i - 1; j >= 0; --j) {
		if(groups[j].second.front()->chrom() != chr)
			return NULL;
		else if(groups[j].first == type)
			return groups[j].second.back();
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_next_record(
								VCFFillableRecord::RecordType type,
								int i, const vector<Group>& groups) const {
	const string&	chr = groups[i].second.front()->chrom();
	for(int j = i + 1; j < (int)groups.size(); ++j) {
		if(groups[j].second.front()->chrom() != chr)
			return NULL;
		else if(groups[j].first == type)
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

// 親のどちらの染色体から来ているかの確率
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

vector<double> VCFFillable::probs_from_which_chrom(RecordSet& rs,
												size_t i, bool is_mat) const {
	if((is_mat && rs.record->is_mat_homo()) ||
					(!is_mat && rs.record->is_pat_homo()))
		return vector<double>(2U, 0.5);
	else if(is_mat)
		return probs_from_which_chrom(rs.prev_mat_from(i), rs.next_mat_from(i));
	else
		return probs_from_which_chrom(rs.prev_pat_from(i), rs.next_pat_from(i));
}

double VCFFillable::likelihood_each(const string& gt,
									const vector<double>& probs_mat,
									const vector<double>& probs_pat,
									int mat_phasing, int pat_phasing) const {
	const int	sum = Genotype::sum_gt(gt);
	double	likelihood = 0.0;
	for(int k = 0; k < 4; ++k) {
		const int	i = k >> 1;		// 母親は0|1か1|0か
		const int	j = k & 1;
		if((((i + mat_phasing) & 1) + ((j + pat_phasing) & 1)) == sum)
			likelihood += probs_mat[i] * probs_pat[j];
	}
	return log(likelihood);
}



double VCFFillable::compute_phasing_likelihood_each(RecordSet& rs, size_t i,
									int mat_phasing, int pat_phasing) const {
	const auto	probs_mat = probs_from_which_chrom(rs, i, true);
	const auto	probs_pat = probs_from_which_chrom(rs, i, false);
	return likelihood_each(rs.gt(i), probs_mat, probs_pat,
									mat_phasing, pat_phasing);
}

// 親のGTがひっくり返っているかひっくり返っていないのか仮定して尤度を計算する
// GTが0/1のとき、phasingが0なら0|1、1なら1|0
// GTが0|0のとき、どちらでも0|0
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
	
	if(rs.record->is_mat_homo()) {
		const string&	gt = rs.record->mat_gt();
		rs.record->set_mat_GT(gt.substr(0, 1) + '|' + gt.substr(2, 1));
	}
	else
		rs.record->set_mat_GT(mat_phasing == 0 ? "0|1" : "1|0");
	if(rs.record->is_pat_homo()) {
		const string&	gt = rs.record->pat_gt();
		rs.record->set_pat_GT(gt.substr(0, 1) + '|' + gt.substr(2, 1));
	}
	else
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
		if(!record->is_homo(0) || !record->is_homo(1)) {
			RecordSet	record_set(record, prev_mat_record, next_mat_record,
											prev_pat_record, next_pat_record);
			determine_phasing(record_set);
			impute_core(record_set);
		}
	}
}

template<typename Iter>
VCFFillableRecord *VCFFillable::find_neighbor_same_type_record(
							size_t i, size_t c, Iter first, Iter last) const {
	const VCFFillableRecord::RecordType	type = fillable_records[i]->get_type();
	const string&	chromosome = fillable_records[i]->chrom();
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
	
	return find_neighbor_same_type_record(i, c, fillable_records.rend() - i + 1,
													fillable_records.rend());
}

VCFFillableRecord *VCFFillable::find_next_same_type_record(
													size_t i, size_t c) const {
	if(i == fillable_records.size() - 1)
		return NULL;
	
	return find_neighbor_same_type_record(i, c,
											fillable_records.begin() + i + 1,
											fillable_records.end());
}

const VCFFillable::RecordSet *VCFFillable::create_recordset(
										size_t i, size_t c, bool is_mat) const {
	auto	*record = fillable_records[i];
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
	
	// patはホモだから1でも2でもよい
	vector<Pair>	pairs;	// [(mat_from, pat_from)]
	if(rs->prev_mat_record != NULL)
		pairs.push_back(Pair(prev_mat_from, 1));
	if(rs->next_mat_record != NULL && next_mat_from != prev_mat_from)
		pairs.push_back(Pair(next_mat_from, 1));
	if(pairs.empty())
		return;
	
	const int	mat_from = select_mat(pairs, rs);
	rs->record->set_GT(c, rs->record->gt_from_parent(mat_from, 1));
	delete rs;
}

void VCFFillable::impute_NA_mat(size_t i) {
	auto	*record = fillable_records[i];
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
	
	vector<Pair>	pairs;
	if(rs->prev_pat_record != NULL)
		pairs.push_back(Pair(1, prev_pat_from));
	if(rs->next_pat_record != NULL && next_pat_from != prev_pat_from)
		pairs.push_back(Pair(1, next_pat_from));
	if(pairs.empty())
		return;
	
	const int	pat_from = select_pat(pairs, rs);
	rs->record->set_GT(c, rs->record->gt_from_parent(1, pat_from));
	delete rs;
}

void VCFFillable::impute_NA_pat(size_t i) {
	auto	*record = fillable_records[i];
	for(size_t c = 11U; c != samples.size() + 9; ++c) {
		if(record->get_GT(c-9) == "./.")
			impute_NA_pat_each(i, c);
	}
}

void VCFFillable::impute() {
cout << get_samples()[0] << " " << get_samples()[1] << endl;
	for(auto p = fillable_records.begin(); p != fillable_records.end(); ++p)
		(*p)->modify();
	
	vector<Group>	groups = group_records();
	for(size_t i = 0U; i < groups.size(); ++i) {
const auto	group = groups[i].second;
for(auto p = group.begin(); p != group.end(); ++p) {
auto record = *p;
if(record->pos() == 183065) {
cout << record->pos() << endl;
}
}
		if(groups[i].second.front()->is_filled_type())
			this->phase(i, groups);
	}
	
	for(size_t i = 0U; i < fillable_records.size(); ++i) {
		auto	*record = fillable_records[i];
		if(record->is_mat_type())
			impute_NA_mat(i);
		else if(record->is_pat_type())
			impute_NA_pat(i);
	}
	
	for(auto p = fillable_records.begin(); p != fillable_records.end(); ++p)
		(*p)->fill_PGT();
}

VCFFillable *VCFFillable::convert(const VCFFamily *vcf) {
	vector<VCFFillableRecord *>	new_records;
	const STRVEC&	samples = vcf->get_samples();
	VCFFillable	*new_vcf = new VCFFillable(vcf->get_header(),
												samples, new_records);
	vcf->copy_chrs(new_vcf);
	for(size_t i = 0U; i < vcf->size(); ++i) {
		const VCFFamilyRecord	*record = vcf->get_record(i);
		auto	new_record = VCFFillableRecord::from_VCFFamilyRecord(record,
																	samples);
		new_records.push_back(new_record);
	}
	new_vcf->set_records(new_records);
	return new_vcf;
}

void VCFFillable::determine_parents_type(VCFFillableRecord *record,
									const TypeDeterminer& determiner) const {
	const tuple<int,int,int> counter = record->count_gts();
	// pは2つのint_gtを4進数でパッキングしたもの
	const int	p = determiner.determine(counter);
	record->determine_parents_type(p);
}

void VCFFillable::replace_filled_records(
					const vector<const VCFRecord *>& orig_records) {
	TypeDeterminer	determiner((int)samples.size(), 0.05);
	for(size_t i = 0U; i < size(); ++i) {
		const VCFRecord	*orig_record = orig_records[i];
		VCFFillableRecord	*record = fillable_records[i];
		if(record->get_v().size() == 9U) {
			auto	v = orig_record->extract_v(record->get_samples());
			record->set(v, record->get_type());
			// 親が間違っているか、ここでチェックする
			// extract_VCFsでチェックするとコストが大きい
if(record->pos() == 183065)
cout << record->pos() << endl;
			determine_parents_type(record, determiner);
		}
	}
}

void VCFFillable::replace_in_thread(void *config) {
	const auto	*c = (ConfigReplaceThread *)config;
	const auto&	vcfs = c->vcfs;
	for(size_t i = c->first; i < vcfs.size(); i += c->num_thread) {
		vcfs[i]->replace_filled_records(c->orig_records);
	}
}

void VCFFillable::replace_filled_records(const vector<VCFFillable *>& vcfs,
													VCFHuge *orig_vcf, int T) {
	// 穴埋めをマルチスレッドで動かすために、必要な分だけ読む
	vector<const VCFRecord *>	orig_records;
	VCFFillable	*vcf = vcfs.front();
	for(size_t i = 0U; i < vcf->size(); ++i) {
		const POSITION	pos = vcf->record_position(*(vcf->get_record(i)));
		const VCFRecord	*record = orig_vcf->proceed(pos);
		orig_records.push_back(record);
	}
	
	vector<ConfigReplaceThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigReplaceThread(vcfs, orig_records, i, T);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i) {
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&replace_in_thread, (void *)configs[i]);
	}
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		replace_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	Common::delete_all(orig_records);
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

void VCFFillable::impute_in_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	const auto&	vcfs = c->vcfs;
	for(size_t i = c->first; i < vcfs.size(); i += c->num_thread) {
		vcfs[i]->impute();
	}
}

void VCFFillable::impute_all_in_multithreads(
							const vector<VCFFillable *>& vcfs, int T) {
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(vcfs, i, T);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i) {
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&impute_in_thread, (void *)configs[i]);
	}
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	for(int i = 0; i < T; ++i)
		delete configs[i];
}

VCFSmall *VCFFillable::join_vcfs(const vector<VCFFamily *>& vcfs_,
											const string& path_VCF, int T) {
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
	VCFFillable::replace_filled_records(filled_vcfs, orig_vcf, T);
#ifndef DEBUG
	VCFFillable::impute_all_in_multithreads(filled_vcfs, T);
#else
	for(auto p = filled_vcfs.begin(); p != filled_vcfs.end(); ++p) {
		(*p)->impute();
	}
#endif
	
	VCFSmall	*merged_vcf = merge_vcfs(filled_vcfs);
	Common::delete_all(filled_vcfs);
	delete orig_vcf;
	return merged_vcf;
}
