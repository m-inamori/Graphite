#include <cmath>
#include "../include/SelfParentImputerLessImputed.h"
#include "../include/common.h"

using namespace std;

SelfParentImputerLessImputed::SelfParentImputerLessImputed(
										const vector<VCFRecord *>& rs,
									 	const vector<vector<int>>& ref_hs,
										const Map& map_, double w) :
					VCFHMM(rs, map_, w), records(rs), ref_haps(ref_hs),
					prev_h_table(collect_possible_previous_hidden_states()),
					Cp(calc_Cp(rs)),
					Epc{{log(1.0-w*2),  log(w/2),     log(w/2),      log(w)},
						{log(0.25-w/4), log(0.5-w/2), log(0.25-w/4), log(w)},
						{log(0.25-w/4), log(0.5-w/2), log(0.25-w/4), log(w)},
						{log(w/2),      log(w/2),     log(1.0-w*2),  log(w)}}
					{ }

vector<double> SelfParentImputerLessImputed::calc_Cp(
							const vector<VCFRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

vector<vector<int>>
SelfParentImputerLessImputed::collect_possible_previous_hidden_states() {
	const size_t	L = num_states();
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hp1 = h % NH();
		const int	hp2 = h / NH();
		// hp2 is the same
		for(int h1 = 0; h1 < (int)NH(); ++h1) {
			const int	prev_h1 = h1 + hp2 * NH();
			prev_h_table[h].push_back(prev_h1);
		}
		// hp1 is the same
		for(int h2 = 0; h2 < (int)NH(); ++h2) {
			if(h2 == hp2)	// for duplication
				continue;
			const int	prev_h1 = hp1 + h2 * NH();
			prev_h_table[h].push_back(prev_h1);
		}
	}
	
	return prev_h_table;
}

double SelfParentImputerLessImputed::progs_emission_probability(
														const vector<int>& ocs,
														int parent_gt) const {
	double	Ec = 0.0;
	for(auto p = ocs.begin(); p != ocs.end(); ++p) {
		Ec += Epc[parent_gt][*p];
	}
	return Ec;
}

double SelfParentImputerLessImputed::emission_probability(
												size_t i, int h, int op,
												const vector<int>& ocs) const {
	const int	gt_parent = parent_genotype(h, i);
	const double	Ep = E[gt_parent][op];	// parent emission
	// progenies emission
	const double	Ec = progs_emission_probability(ocs, gt_parent);
	return Ep + Ec;
}

double SelfParentImputerLessImputed::transition_probability(size_t i,
													int prev_h, int h) const {
	const int	prev_hp1 = prev_h % NH();
	const int	prev_hp2 = prev_h / NH();
	const int	hp1 = h % NH();
	const int	hp2 = h / NH();
	const double	cp = Cp[i-1];
	return (log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
			log(prev_hp2 != hp2 ? cp : 1.0 - cp));
}

vector<SelfParentImputerLessImputed::DP>
						SelfParentImputerLessImputed::initialize_dp() const {
	const size_t	L = num_states();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFRecord	*record = records[0];
	const int	op = Genotype::gt_to_int(record->get_gt(0));
	vector<int>	ocs;	// observed progenies' genotypes
	for(size_t ic = 0; ic != num_progenies(); ++ic) {
		ocs.push_back(record->get_int_gt(ic+1));
	}
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(0, h, op, ocs);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void SelfParentImputerLessImputed::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = num_states();
	const VCFRecord	*record = records[i];
	const int	op = Genotype::gt_to_int(record->get_gt(0));
	vector<int>	ocs;	// observed progenies' genotypes
	for(size_t ic = 0; ic != num_progenies(); ++ic) {
		ocs.push_back(Genotype::get_int_gt(record->get_gt(ic+1)));
	}
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, op, ocs);
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void SelfParentImputerLessImputed::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		const int	parent_gt = parent_genotype(hs[i], i);
		records[i]->set_GT(0, Genotype::int_to_phased_gt(parent_gt));
	}
}

void SelfParentImputerLessImputed::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
