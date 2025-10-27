#include <sstream>
#include "../include/VCFHeteroImpHomo.h"

using namespace std;


//////////////////// VCFHeteroImpHomo ////////////////////

VCFHeteroImpHomo::DP VCFHeteroImpHomo::init_dp() const {
	const int	n = num_progenies();
	const int	L = (n*2+1) << n;
	DP	dp(L, Value::create_default(n));
	dp[0] = Value(0, State(0, n), 0);
	return dp;
}

// When heterozygous parents are determined
// are the genotype of each progeny correct even if it does not crossover?
vector<bool> VCFHeteroImpHomo::is_right_gt(int order, int gt_imputed,
							const State& state, const GenoRecord *record) const {
	const size_t	n = num_progenies();
	vector<bool>	bs(n, false);
	for(size_t i = 0; i < n; ++i) {
		const int	gt = record->unphased(i+2);
		if(Genotype::is_NA(gt)) {
			// If it is N/A, do not crossover
			// considering that there is no mistake
			bs[i] = true;
		}
		else {
			const int	prev_h = state.haplotype(i);
			const int	new_gt = gt_by_haplotype(prev_h, gt_imputed, order);
			bs[i] = new_gt == gt;
		}
	}
	return bs;
}

vector<pair<VCFHeteroImpHomo::State, VCFHeteroImpHomo::Value>>
VCFHeteroImpHomo::next_states(State state, const GenoRecord *record) const {
	vector<pair<State, Value>>	states;
	const int	n = (int)num_progenies();
	const int	gt_imputed = record->unphased(imputed_index());
	for(int order = 0; order < 2; ++order) {
		const vector<bool>	bs = is_right_gt(order, gt_imputed, state, record);
		int	num_falses = 0;
		for(auto p = bs.begin(); p != bs.end(); ++p) {
			if(!*p)
				num_falses += 1;
		}
		State	s0 = State(0, n);
		s0.set_num_crossovers(state.num_crossovers());
		for(int s1 = 0; s1 < (1 << num_falses); ++s1) {
			State	s(s0.s, n);
			int	mis = 0;
			int	j = 0;	// Increment until it reaches num_falses
			bool	over_crossovers = false;
			for(int i = 0; i < n; ++i) {
				int	h;
				const int	prev_h = state.haplotype(i);
				if(!bs[i]) {	// if no crossover, wrong genotype
					h = (s1 >> j) & 1;
					if(h != prev_h) {	// crossover
						if(s.is_full_crossovers()) {
							over_crossovers = true;
							break;
						}
						s.increment_num_crossovers();
					}
					else {
						mis += 1;
					}
					j += 1;
				}
				else {
					h = prev_h;
				}
				s.set_haplotype(h, i);
			}
			if(!over_crossovers)
				states.push_back(make_pair(s, Value(mis, state, order)));
		}
	}
	return states;
}

VCFHeteroImpHomo::DP VCFHeteroImpHomo::update_dp(const DP& dp,
											const GenoRecord *record) const {
	const int	n = num_progenies();
	const int	L = (n*2+1) << n;
	DP	new_dp(L, Value(INF, State(0, n), 0));
	for(int s = 0; s < L; ++s) {
		if(!dp[s].is_valid())
			continue;
		State	state(s, n);
		const auto	v = next_states(state, record);
		for(auto p = v.begin(); p != v.end(); ++p) {
			const State	new_state = p->first;
			const Value	new_value = p->second;
			const Value	value = dp[s].add(new_value, state);
			new_dp[new_state.s] = min(new_dp[new_state.s], value);
		}
	}
	return new_dp;
}

void VCFHeteroImpHomo::trace_back(State state, const vector<DP>& dps) {
	for(size_t k = records.size(); k > 0; --k) {
		const size_t	i = k - 1;
		const int	order = dps[i+1][state.s].order;
		GenoRecord	*record = records[i];
		const int	hetero_GT = order == 0 ? Genotype::PH_01 : Genotype::PH_10;
		record->set_geno(non_imputed_index(), hetero_GT);
		const char	gt_imp = record->get_allele(imputed_index(), 0);
		for(size_t j = 2; j < get_samples().size(); ++j) {
			const int	h = state.haplotype(j-2);
//if(j >= record->get_geno().size())
//cerr << "j : " << j << " size : " << get_samples().size() << " " << record->get_geno().size() << endl;
			if(is_mat_hetero)
				record->set_geno(j, Genotype::from_alleles(order^h, gt_imp));
			else
				record->set_geno(j, Genotype::from_alleles(gt_imp, order^h));
		}
		state = dps[k][state.s].state;
	}
}

void VCFHeteroImpHomo::impute() {
	vector<DP>	dps(1, init_dp());
	for(auto p = records.begin(); p != records.end(); ++p)
		dps.push_back(update_dp(dps.back(), *p));
	
	pair<int, int>	min_pair = make_pair(dps.back()[0].num_mis, 0);
	for(int state = 1; state < (int)dps.back().size(); ++state) {
		const pair<int, int>	p(dps.back()[state].num_mis, state);
		if(p < min_pair)
			min_pair = p;
	}
	const int	min_state = min_pair.second;
	trace_back(State(min_state, 0), dps);
}
