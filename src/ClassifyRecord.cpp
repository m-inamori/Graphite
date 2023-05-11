#include <numeric>
#include <algorithm>
#include "../include/ClassifyRecord.h"
#include "../include/VCFFamily.h"

using namespace std;


//////////////////// ClassifyRecord ////////////////////

ClassifyRecord *ClassifyRecord::instance = nullptr;

ClassifyRecord::~ClassifyRecord() {
	for(auto p = memo.begin(); p != memo.end(); ++p) {
		delete p->second;
	}
}

ClassifyRecord *ClassifyRecord::get_instance() {
	if(!instance)
		instance = new ClassifyRecord();
	return instance;
}

pair<ParentComb, WrongType> ClassifyRecord::classify(
									const VCFFamilyRecord *record,
									const TypeDeterminer *td) {
	const vector<int>	gts = record->get_progeny_int_gts();
	const tuple<int, int, int>	counter = count_int_gts(gts);
	const auto	combs = td->determine(counter);
	return classify_record_core(combs, record->mat_int_gt(),
												record->pat_int_gt());
}

const TypeDeterminer *ClassifyRecord::get_TypeDeterminer(size_t n, double a) {
	auto	p = memo.find(make_pair(n, a));
	if(p != memo.end())
		return p->second;
	
	const TypeDeterminer	*td = new TypeDeterminer(n, a);
	memo.insert(make_pair(make_pair(n, a), td));
	return td;
}

tuple<size_t, size_t, size_t> ClassifyRecord::count_int_gts(
											const vector<int>& gts) const {
	vector<size_t>	ns(3, 0);
	for(auto p = gts.begin(); p != gts.end(); ++p) {
		if(*p != -1)
			ns[*p] += 1;
	}
	return tuple<size_t, size_t, size_t>(ns[0], ns[1], ns[2]);
}

vector<ClassifyRecord::GTComb> ClassifyRecord::filter_pairs(
											const vector<GTComb>& combs) const {
	if(combs.size() == 1)
		return combs;
	
	// 親の組合せの候補が複数ある場合、取れる確率が大きく違えば一つにする
	// p0(1-p1)(1-p2)などを確率のようなものとみなす
	vector<double>	ps(combs.size());	// p0(1-p1)(1-p2)などを格納する
	for(size_t i = 0; i < combs.size(); ++i) {
		double	p = 1.0;
		for(size_t k = 0; k < combs.size(); ++k)
			p *= i == k ? 1.0 - combs[k].first : combs[k].first + 0.01;
		ps[i] = p;
	}
	
	double s = std::accumulate(ps.begin(), ps.end(), 0.0);
	for(size_t i = 0; i < ps.size(); ++i) {
		if(ps[i] / s >= 0.99)
			return vector<GTComb>(1, combs[i]);
	}
	return combs;
}

WrongType ClassifyRecord::select_wrong_type(ParentComb comb,
											int mat_gt, int pat_gt) const {
	if(TypeDeterminer::is_same_parent_gts(comb)) {
		const int	gt = TypeDeterminer::int_gt_pair(comb).first;
		if(mat_gt == gt && pat_gt == gt)
			return WrongType::RIGHT;
		else
			return WrongType::MODIFIABLE;
	}
	else {
		// 0/0 x 0/1 -> 2, 0/0 x 1/1 -> 1, 0/1 x 1/1 -> 0
		const int	avoiding_gt = TypeDeterminer::get_avoiding_gt(comb);
		if(mat_gt == pat_gt)
			return WrongType::UNMODIFIABLE;
		else if(mat_gt == -1 && (pat_gt != -1 && pat_gt != avoiding_gt))
			return WrongType::MODIFIABLE;
		else if(pat_gt == -1 && (mat_gt != -1 && mat_gt != avoiding_gt))
			return WrongType::MODIFIABLE;
		else if(mat_gt != avoiding_gt && pat_gt != avoiding_gt)
			return WrongType::RIGHT;
		else
			return WrongType::MODIFIABLE;
	}
}

bool ClassifyRecord::is_matched(int mat_gt, int pat_gt, ParentComb comb) const {
	const auto	gt_pair = TypeDeterminer::int_gt_pair(comb);
	return (mat_gt == gt_pair.first && pat_gt == gt_pair.second) ||
			(pat_gt == gt_pair.first && mat_gt == gt_pair.second);
}

std::pair<ParentComb, WrongType> ClassifyRecord::select_pair(
												vector<GTComb>& combs,
												int mat_gt, int pat_gt) const {
	// combsの中に確率が大きいものがあり、それ以外は小さければ、それだけにする
	if(combs.size() >= 2) {
		std::sort(combs.begin(), combs.end());	// 確率でソート
		vector<double>	ps;
		for(auto p = combs.begin(); p != combs.end(); ++p) {
			ps.push_back(p->first);
		}
		vector<double>	qs;
		for(auto p = ps.begin(); p != ps.end(); ++p) {
			double	prob = 1.0;
			for(auto q = ps.begin(); q != ps.end(); ++q)
				prob *= p == q ? 1.0 - *q : *q + 0.01;
			qs.push_back(prob);
		}
		if(qs.front() / std::accumulate(qs.begin(), qs.end(), 0.0) >= 0.99)
			combs = vector<GTComb>(1, combs.front());
	}
	
	const ParentComb	comb = combs.front().second;
	if(is_matched(mat_gt, pat_gt, comb))
		return make_pair(comb, WrongType::RIGHT);
	else if(mat_gt == pat_gt)	// どちらが正しいのかわからない
		return make_pair(ParentComb::PNA, WrongType::MIX);
	
	// 最も優先順位が高いペアだけ片方だけ合っているなら修正可能
	vector<bool>	bs;
	for(auto p = combs.begin(); p != combs.end(); ++p) {
		const auto	gt_pair = TypeDeterminer::int_gt_pair(p->second);
		const bool	b = (mat_gt == gt_pair.first || mat_gt == gt_pair.second ||
						 pat_gt == gt_pair.first || pat_gt == gt_pair.second);
		if((p == combs.begin()) ^ b)
			return make_pair(ParentComb::PNA, WrongType::MIX);
	}
	return make_pair(comb, WrongType::MODIFIABLE);
}

pair<ParentComb, WrongType> ClassifyRecord::classify_record_core(
												const vector<GTComb>& pairs_,
												int mat_gt, int pat_gt) const {
	if(pairs_.size() == 0)
		return make_pair(ParentComb::PNA, WrongType::UNSPECIFIED);
	
	auto	combs = filter_pairs(pairs_);
	if(combs.size() == 1) {
		const ParentComb	p = combs.front().second;
		return make_pair(p, select_wrong_type(p, mat_gt, pat_gt));
	}
	else {
		return select_pair(combs, mat_gt, pat_gt);
	}
}
