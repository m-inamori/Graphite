#include <cmath>
#include "../include/SelfParentImputer.h"
#include "../include/common.h"

using namespace std;

SelfParentImputer::SelfParentImputer(const vector<GenoRecord *>& rs,
									 const vector<vector<int>>& ref_hs,
									 size_t iprog, const Map& map_, double w) :
			VCFHMM(rs, map_, w), ref_records(rs), ref_haps(ref_hs),
			ic(iprog), prev_h_table(collect_possible_previous_hidden_states()),
			Cc(calc_Cc(rs)), Cp(calc_Cp(rs)),
			Epc{{log(1.0-w*2),  log(w/2),     log(w/2),      log(w)},
				{log(0.25-w/4), log(0.5-w/2), log(0.25-w/4), log(w)},
				{log(0.25-w/4), log(0.5-w/2), log(0.25-w/4), log(w)},
				{log(w/2),      log(w/2),     log(1.0-w*2),  log(w)}}
			{ }

vector<double> SelfParentImputer::calc_Cc(
							const vector<GenoRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> SelfParentImputer::calc_Cp(
							const vector<GenoRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

tuple<int, int, int, int> SelfParentImputer::decode_state(int h) const {
	// 1st bit: Which parent's haplotype does the left side of the progeny come from?
	// 2nd bit: Which parent's haplotype does the right side of the progeny come from?
	// At least one of the parent's haplotypes must be the same as the progeny's
	const int	hc = h & 3;
	const int	hp = h >> 2;
	if(hp < (int)NH() * 2 + 4) {
		const int	hp1 = hp % (NH() + 2);
		const int	hp2 = hp / (NH() + 2);
		return make_tuple(hp1, hp2, hc & 1, hc >> 1);
	}
	else {
		const int	hp1 = (hp - NH() * 2) & 1;
		const int	hp2 = (hp - NH() * 2) >> 1;
		return make_tuple(hp1, hp2, hc & 1, hc >> 1);
	}
}

double SelfParentImputer::phased_emission_probability(int h, size_t i,
															int prog_gt) const {
	const int	phased_prog_gt = progeny_genotype(h, i, prog_gt);
	const int	d = phased_prog_gt ^ prog_gt;
	return log((((d & 1) == 0) ? 0.99 : 0.01) *
				(((d >> 1) == 0) ? 0.99 : 0.01));
}

double SelfParentImputer::non_phased_emission_probability(
												const vector<int>& ocs,
												int parent_gt) const {
	double	s = 0.0;
	for(auto p = ocs.begin(); p != ocs.end(); ++p) {
		s += Epc[parent_gt][*p];
	}
	return s;
}

double SelfParentImputer::emission_probability(int h, size_t i, int op,
												int prog_gt,
												const vector<int>& ocs) const {
	const int	parent_gt = parent_genotype(h, i, prog_gt);
	const double	Ep = E[parent_gt][op];
	const double	Ec = phased_emission_probability(h, i, prog_gt);
	const double	Ecs = non_phased_emission_probability(ocs, parent_gt);
	return Ep + Ec + Ecs;
}

double SelfParentImputer::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cp = Cp[i-1];
	const double	cc = Cc[i-1];
	const auto	t1 = decode_state(h);
	const auto	t2 = decode_state(prev_h);
	return (log(get<0>(t1) != get<0>(t2) ? cp : 1.0 - cp) +
			log(get<1>(t1) != get<1>(t2) ? cp : 1.0 - cp) +
			log(get<2>(t1) != get<2>(t2) ? cc : 1.0 - cc) +
			log(get<3>(t1) != get<3>(t2) ? cc : 1.0 - cc));
}

vector<SelfParentImputer::DP> SelfParentImputer::initialize_dp() const {
	const size_t	L = num_states();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const GenoRecord	*record = records[0];
	const int	op = record->unphased(0);
	const int	prog_gt = record->get_geno()[ic+1] & 3;
	vector<int>	ocs;
	for(size_t j = 0; j < num_progenies(); ++j) {
		if(j != ic) {
			ocs.push_back(record->unphased(j+1));
		}
	}
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(h, 0, op, prog_gt, ocs);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
SelfParentImputer::collect_possible_previous_hidden_states() const {
	const size_t	L = num_states();
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const auto	t1 = decode_state(h);
		for(int prev_h = 0; prev_h < (int)L; ++prev_h) {
			const auto	t2 = decode_state(prev_h);
			int	diff_counter = (get<0>(t1) != get<0>(t2) ? 1 : 0) +
							   (get<1>(t1) != get<1>(t2) ? 1 : 0) +
							   (get<2>(t1) != get<2>(t2) ? 1 : 0) +
							   (get<3>(t1) != get<3>(t2) ? 1 : 0);
			if(diff_counter <= 1) {
				prev_h_table[h].push_back(prev_h);
			}
		}
	}
	
	return prev_h_table;
}

void SelfParentImputer::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = num_states();
	const GenoRecord	*record = records[i];
	const int	op = record->unphased(0);
	const int	prog_gt = record->get_geno()[ic+1] & 3;
	vector<int>	ocs;
	for(size_t j = 0; j < num_progenies(); ++j) {
		if(j != ic) {
			ocs.push_back(record->unphased(j+1));
		}
	}
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(h, i, op, prog_gt, ocs);
		
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void SelfParentImputer::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		GenoRecord	*record = records[i];
		const int	prog_gt = record->get_geno()[ic+1];
		const int	phased_gt = parent_genotype(hs[i], i, prog_gt);
		record->set_geno(0, phased_gt | 4);
	}
}

void SelfParentImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
