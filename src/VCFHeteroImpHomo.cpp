#include <sstream>
#include "../include/VCFHeteroImpHomo.h"

using namespace std;


//////////////////// VCFHeteroImpHomo ////////////////////

VCFHeteroImpHomo::DP VCFHeteroImpHomo::init_dp() const {
	const int	n = num_progenies();
	const int	L = (n*2+1) << (n+1);
	DP	dp(L, DPValue(INF, 0));
	dp[0] = DPValue(0, 0);
	return dp;
}

vector<pair<int, int>> VCFHeteroImpHomo::next_states_small(int state,
												const VCFRecord *record) const {
	vector<pair<int, int>>	states;
	const int	n = num_progenies();
	const int	prev_num_crossovers = get_crossover(state);
	for(int s = 0; s < (1 << (n+1)); ++s) {
		const int	order = get_order(s);
		int	mis = 0;
		int	num_crossovers = prev_num_crossovers;
		const int	gt_imputed = record->get_int_gt(imputed_index());
		for(size_t i = 2; i < (size_t)n+2; ++i) {
			const int	gt = record->get_int_gt(i);
			const int	prev_h = get_hap(state, i-2);
			const int	h = get_hap(s, i-2);
			if(prev_h != h)
				num_crossovers += 1;
			if(gt != -1 && gt != gt_imputed/2 + (order ^ h))
				mis += 1;
		}
		if(num_crossovers <= n * 2) {
			const int	s1 = s | (num_crossovers << (n+1));
			states.push_back(make_pair(s1, mis));
		}
	}
	return states;
}

// When heterozygous parents are determined
// are the genotype of each progeny correct even if it does not crossover?
vector<bool> VCFHeteroImpHomo::is_right_gt(int order, int gt_imputed,
									int state, const VCFRecord *record) const {
	const size_t	n = num_progenies();
	vector<bool>	bs(n, false);
	for(size_t i = 2; i < n+2; ++i) {
		const int	gt = record->get_int_gt(i);
		if(gt == -1) {
			// If it is N/A, do not crossover
			// considering that there is no mistake
			bs[i-2] = true;
		}
		else {
			const int	prev_h = get_hap(state, i-2);
			const int	new_gt = gt_by_haplotype(prev_h, gt_imputed, order);
			bs[i-2] = new_gt == gt;
		}
	}
	return bs;
}

vector<pair<int, int>> VCFHeteroImpHomo::next_states_large(int state,
												const VCFRecord *record) const {
	vector<pair<int, int>>	states;
	const int	n = (int)num_progenies();
	const int	prev_num_crossovers = get_crossover(state);
	const int	gt_imputed = record->get_int_gt(imputed_index());
	for(int order = 0; order < 2; ++order) {
		const vector<bool>	bs = is_right_gt(order, gt_imputed,
													state, record);
		int	num_falses = 0;
		for(auto p = bs.begin(); p != bs.end(); ++p) {
			if(!*p)
				num_falses += 1;
		}
		int	s = order | (prev_num_crossovers << (n+1));
		int	mis = 0;
		for(int s1 = 0; s1 < (1 << num_falses); ++s1) {
			int	j = 0;
			for(int i = 2; i < n+2; ++i) {
				const int	prev_h = get_hap(state, i-2);
				if(!bs[i-2]) {
					const int	h = (s1 >> j) & 1;
					if(h != prev_h) {	// crossover
						if(s >> (n+1) < n * 2) {
							s |= h << (i-1);
							s += 1 << (n+1);
						}
					}
					else {
						mis += 1;
					}
				}
				else {
					s |= prev_h << (i-1);
				}
			}
		}
	}
	return states;
}

vector<pair<int, int>> VCFHeteroImpHomo::next_states(int state,
												const VCFRecord *record) const {
	if(num_progenies() < 5)
		return next_states_small(state, record);
	else
		return next_states_large(state, record);
}

VCFHeteroImpHomo::DP VCFHeteroImpHomo::update_dp(const DP& dp,
												const VCFRecord *record) const {
	const int	n = num_progenies();
	const int	L = (n*2+1) << (n+1);
	DP	new_dp(L, DPValue(INF, 0));
	for(int state = 0; state < L; ++state) {
		const int	mis = dp[state].first;
		const auto	v = next_states(state, record);
		for(auto p = v.begin(); p != v.end(); ++p) {
			const int	new_state = p->first;
			const int	new_mis = p->second;
			const pair<int, int>	value(mis+new_mis, state);
			new_dp[new_state] = min(new_dp[new_state], value);
		}
	}
	return new_dp;
}

void VCFHeteroImpHomo::trace_back(int state, const vector<DP>& dps) {
	for(size_t k = records.size(); k > 0; --k) {
		const size_t	i = k - 1;
		const int	order = get_order(state);
		VCFRecord	*record = records[i];
		string	hetero_GT = order == 0 ? "0|1" : "1|0";
		record->set_GT(non_imputed_index(), hetero_GT);
		const char	gt_imp = record->get_gt(imputed_index()).c_str()[0];
		for(size_t j = 2; j < samples.size(); ++j) {
			const int	h = get_hap(state, j-2);
			stringstream	ss;
			if(is_mat_hetero)
				ss << (order ^ h) << '|' << gt_imp;
			else
				ss << gt_imp << '|' << (order ^ h);
			const string	GT = ss.str();
			record->set_GT(j, GT);
		}
		state = dps[k][state].second;
	}
}

void VCFHeteroImpHomo::impute() {
	vector<DP>	dps(1, init_dp());
	for(auto p = records.begin(); p != records.end(); ++p)
		dps.push_back(update_dp(dps.back(), *p));
	
	int	min_state = 0;
	DPValue	min_value = dps.back()[0];
	for(int state = 1; state < (int)dps.back().size(); ++state) {
		const DPValue	value = dps.back()[state];
		if(value < min_value) {
			min_state = state;
			min_value = value;
		}
	}
	trace_back(min_state, dps);
}
