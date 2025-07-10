#include <cmath>
#include "../include/SelfImputer.h"
#include "../include/common.h"

using namespace std;

SelfImputer::SelfImputer(const vector<VCFRecord *>& rs,
							const vector<vector<int>>& ref_hs,
							const Map& map_, double w) :
					VCFHMM(rs, map_, w), records(rs),
					ref_haps(ref_hs), Cc(calc_Cc(rs))
					{ }

vector<double> SelfImputer::calc_Cc(
							const vector<VCFRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

bool SelfImputer::is_only_one_or_zero_crossover(int hp1, int hp2,
												int hc, int prev_h) const {
	const auto	t = decode_state(prev_h);
	const int	prev_hp1 = get<0>(t);
	const int	prev_hp2 = get<1>(t);
	const int	prev_hc  = get<2>(t);
	if(hp1 != prev_hp1 && hp2 != prev_hp2) {
		return false;
	}
	else if(hp1 == prev_hp1 && hp2 == prev_hp2) {
		int	d = hc ^ prev_hc;
		int	diff_counter = 0;
		while(d > 0) {
			if((d & 1) == 1) {
				diff_counter += 1;
				if(diff_counter > 1)
					return false;
			}
			d >>= 1;
		}
		return true;
	}
	else {	// 片側が乗り換える
		return hc == prev_hc;
	}
}

vector<vector<int>>
SelfImputer::collect_possible_previous_hidden_states() const {
	const size_t	L = num_states();
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {
		const auto	t = decode_state(h);
		const int	hp1 = get<0>(t);
		const int	hp2 = get<1>(t);
		const int	hc  = get<2>(t);
		for(int prev_h = 0; prev_h < (int)L; ++prev_h) {
			if(is_only_one_or_zero_crossover(hp1, hp2, hc, prev_h))
				prev_h_table[h].push_back(prev_h);
		}
	}
	
	return prev_h_table;
}

tuple<int, int, int> SelfImputer::decode_state(int h) const {
	const int	hc = h & ((1 << num_states()) - 1);
	const int	hp = h >> num_states();
	const int	hp1 = hp % NH();
	const int	hp2 = hp / NH();
	return make_tuple(hp1, hp2, hc);
}

vector<int> SelfImputer::compute_progeny_phased_gts(int hc,
													int parent_gt) const {
	vector<int>	prog_gts(num_progenies());
	for(int j = 0; j < (int)num_progenies(); ++j) {
		prog_gts[j] = ((parent_gt >> (hc >> (j*2))) & 1) |
					  ((parent_gt >> (hc >> (j*2+1))) & 1);
	}
	return prog_gts;
}

double SelfImputer::progs_emission_probability(int hc, const vector<int>& ocs,
														int parent_gt) const {
	const auto	phased_gts = compute_progeny_phased_gts(hc, parent_gt);
	double	Ec = 0.0;
	for(size_t j = 0; j < num_progenies(); ++j) {
		Ec += E[phased_gts[j]][ocs[j]];
	}
	return Ec;
}

double SelfImputer::emission_probability(int i, int h, int op,
											std::vector<int>& ocs) const {
	const auto	t = decode_state(h);
	const int	hp1 = get<0>(t);
	const int	hp2 = get<1>(t);
	const int	hc  = get<2>(t);
	const int	gt_parent = parent_genotype(hp1, hp2, i);
	const double	Ep = E[gt_parent][op];	// parent emission
	// progenies emission
	const double	Ec = progs_emission_probability(hc, ocs, gt_parent);
	return Ep + Ec;
}

double SelfImputer::parent_transition_probability(int i, int prev_hp1,
													int prev_hp2,
													int hp1, int hp2) const {
	const double	cp = Cp[i-1];
	return log(hp1 != prev_hp1 ? cp : 1.0 - cp) +
		   log(hp2 != prev_hp2 ? cp : 1.0 - cp);
}

double SelfImputer::progeny_transition_probability(int i, int prev_hc,
																int hc) const {
	const double	cc = Cc[i-1];
	double	T = 0.0;
	for(int j = 0; j < (int)num_progenies()*2; ++j) {
		if(((hc >> j) & 1) != ((prev_hc >> j) & 1))
			T += log(cc);
		else
			T += log(1.0 - cc);
	}
	return T;
}

double SelfImputer::transition_probability(size_t i, int prev_h, int h) const {
	const auto	t1 = decode_state(prev_h);
	const int	prev_hp1 = get<0>(t1);
	const int	prev_hp2 = get<1>(t1);
	const int	prev_hc  = get<2>(t1);
	const auto	t2 = decode_state(h);
	const int	hp1 = get<0>(t2);
	const int	hp2 = get<1>(t2);
	const int	hc  = get<2>(t2);
	const double	Tp = parent_transition_probability(i, prev_hp1, prev_hp2,
																	hp1, hp2);
	const double	Tc = progeny_transition_probability(i, prev_hc, hc);
	return Tp + Tc;
}

vector<SelfImputer::DP> SelfImputer::initialize_dp() const {
	const size_t	L = num_states();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFRecord	*record = records[0];
	const int	op = Genotype::gt_to_int(record->get_gt(0));
	// observed progs
	vector<int>	ocs;
	for(size_t i = 0; i != num_progenies(); ++i) {
		ocs.push_back(Genotype::gt_to_int(record->get_gt(i+1)));
	}
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(0, h, op, ocs);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void SelfImputer::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = num_states();
	const VCFRecord	*record = records[i];
	const int	op = Genotype::gt_to_int(record->get_gt(0));
	// observed progs
	vector<int>	ocs;
	for(size_t i = 0; i != num_progenies(); ++i) {
		ocs.push_back(Genotype::gt_to_int(record->get_gt(i+1)));
	}
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, op, ocs);
		
		for(int prev_h = 0; prev_h < (int)L; ++prev_h) {
			const double	T_all = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (T_all + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void SelfImputer::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		VCFRecord	*record = records[i];
		const auto	t = decode_state(hs[i]);
		const int	hp1 = get<0>(t);
		const int	hp2 = get<1>(t);
		const int	hc  = get<2>(t);
		const int	parent_gt = parent_genotype(hp1, hp2, i);
		record->set_GT(0, Genotype::int_to_phased_gt(parent_gt));
		const auto	prog_gts = compute_progeny_phased_gts(hc, parent_gt);
		for(size_t j = 0; j < num_progenies(); ++j) {
			const string	GT = Genotype::int_to_phased_gt(prog_gts[j]);
			record->set_GT(j + 1, GT);
		}
	}
}

void SelfImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
