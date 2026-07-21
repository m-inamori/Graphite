#include <array>
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
									const TypeDeterminer *td,
									bool one_parent) {
	const vector<int>	gts = record->unphased_gts();
	const array<int, 4>	counter = count_int_gts(gts.begin() + 2, gts.end());
	const auto	combs = td->determine(gts[0], gts[1], counter);
	return classify_record_core(combs, record->unphased_mat(),
										record->unphased_pat(), one_parent);
}

pair<ParentComb, WrongType> ClassifyRecord::classify_self_record(
												const GenoRecord *record,
												const TypeDeterminer *td) {
	const vector<int>	gts = record->unphased_gts();
	const array<int, 4>	counter = count_int_gts(gts.begin() + 1, gts.end());
	const auto	pairs = td->determine(gts[0], gts[0], counter);
	if(pairs.empty())
		return make_pair(ParentComb::PNA, WrongType::UNSPECIFIED);
	
	const int	parent_gt = gts[0];
	const ParentComb	comb = pairs[0].first;
	if(comb == ParentComb::P01x01 &&
						(parent_gt == 1 || Genotype::is_NA(parent_gt)))
		return make_pair(ParentComb::P01x01, WrongType::RIGHT);
	else if(comb == ParentComb::P00x00 &&
						(parent_gt == 0 || Genotype::is_NA(parent_gt)))
		return make_pair(ParentComb::P00x00, WrongType::RIGHT);
	else if(comb == ParentComb::P11x11 &&
						(parent_gt == 2 || Genotype::is_NA(parent_gt)))
		return make_pair(ParentComb::P11x11, WrongType::RIGHT);
	else
		return make_pair(ParentComb::PNA, WrongType::UNSPECIFIED);
}

const TypeDeterminer *ClassifyRecord::get_TypeDeterminer(size_t n) {
	auto	p = memo.find(n);
	if(p != memo.end())
		return p->second;
	
	const TypeDeterminer	*td = new TypeDeterminer(n);
	memo.insert(make_pair(n, td));
	return td;
}

// count the numbers of 0/0, 0/1, 1/1, and ./.
array<int, 4> ClassifyRecord::count_int_gts(const vector<int>& gts) const {
	return count_int_gts(gts.begin(), gts.end());
}

array<int, 4> ClassifyRecord::count_int_gts(
									vector<int>::const_iterator first,
									vector<int>::const_iterator last) const {
	array<int, 4>	ns = {0, 0, 0, 0};
	for(auto p = first; p != last; ++p) {
		ns[*p] += 1;
	}
	return ns;
}

WrongType ClassifyRecord::select_wrong_type(ParentComb comb, int mat_gt,
											int pat_gt, bool one_parent) const {
	if(TypeDeterminer::is_same_parent_genotype(comb)) {
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
		else if((Genotype::is_NA(mat_gt) &&
						(!Genotype::is_NA(pat_gt) && pat_gt != avoiding_gt)) ||
				(Genotype::is_NA(pat_gt) &&
						(!Genotype::is_NA(mat_gt) && mat_gt != avoiding_gt))) {
			return one_parent ? WrongType::RIGHT : WrongType::MODIFIABLE;
		}
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

pair<ParentComb, WrongType> ClassifyRecord::classify_record_core(
												const vector<GTComb>& pairs,
												int mat_gt, int pat_gt,
												bool one_parent) const {
	if(pairs.size() == 1) {
		const ParentComb	p = pairs.front().first;
		return make_pair(p, select_wrong_type(p, mat_gt, pat_gt, one_parent));
	}
	else {
		return make_pair(ParentComb::PNA, WrongType::UNSPECIFIED);
	}
}

pair<ParentComb, FillType> ClassifyRecord::classify_family_record(
												const VCFFamilyRecord *record) {
	if(record->is_NA(0) || record->is_NA(1))
		return make_pair(ParentComb::PNA, FillType::IMPUTABLE);
	
	if(record->is_mat_hetero()) {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P01x01, FillType::IMPUTABLE);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x01, FillType::MAT);
		else
			return make_pair(ParentComb::P01x11, FillType::MAT);
	}
	else if(record->is_00(0)) {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P00x01, FillType::PAT);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x00, FillType::FILLED);
		else
			return make_pair(ParentComb::P00x11, FillType::FILLED);
	}
	else {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P01x11, FillType::PAT);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x11, FillType::FILLED);
		else
			return make_pair(ParentComb::P11x11, FillType::FILLED);
	}
}
