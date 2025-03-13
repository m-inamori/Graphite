#include <cmath>
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

OrphanImputer::OrphanImputer(const vector<VCFRecord *>& rs,
							 const vector<vector<int>>& ref_hs,
							 const Map& map_, double w) :
				VCFHMM(rs, map_, w), ref_records(rs), ref_haps(ref_hs),
				prev_h_table(collect_possible_previous_hidden_states()),
				Cp(calc_Cp(rs)) { }

vector<double> OrphanImputer::calc_Cp(const vector<VCFRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

double OrphanImputer::emission_probability(size_t i, int h,
												int orphan_gt) const {
	const int	phased_gt = compute_phased_gt_by_refhaps(h, i);
	return E[phased_gt][orphan_gt];
}

vector<OrphanImputer::DP> OrphanImputer::initialize_dp(size_t io) const {
	const size_t	L = NH() * NH();
	vector<DP>	dp(M(), DP(L, make_pair(MIN_PROB, 0)));
	const auto	*record = ref_records[0];
	const int	orphan_gt = Genotype::gt_to_int(record->get_gt(io));
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, h, orphan_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
OrphanImputer::collect_possible_previous_hidden_states() const {
	const size_t	L = NH() * NH();
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hp1 = h % NH();
		const int	hp2 = h / NH();
		// Neither side will crossover
		for(int h1 = 0; h1 < (int)NH(); ++h1) {
			const int	prev_h1 = h1 + hp2 * NH();
			prev_h_table[h].push_back(prev_h1);
		}
		for(int h2 = 0; h2 < (int)NH(); ++h2) {
			if(h2 == hp2)	// for duplication
				continue;
			const int	prev_h1 = hp1 + h2 * NH();
			prev_h_table[h].push_back(prev_h1);
		}
	}
	
	return prev_h_table;
}

void OrphanImputer::update_dp(size_t i, size_t io, vector<DP>& dp) const {
	const size_t	L = NH() * NH();
	const VCFRecord	*record = ref_records[i];
	
	// observed
	const int	orphan_gt = Genotype::gt_to_int(record->get_gt(io));
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, orphan_gt);
		
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const double	To = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (To + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void OrphanImputer::update_genotypes(const vector<int>& hs, size_t io) {
	for(size_t i = 0; i < this->M(); ++i) {
		const int	phased_gt = compute_phased_gt_by_refhaps(hs[i], i);
		ref_records[i]->set_GT(io, Genotype::int_to_phased_gt(phased_gt));
	}
}

void OrphanImputer::impute(size_t io) {
	// DP
	vector<DP>	dp = initialize_dp(io);
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, io, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs, io);
}
