#include <cmath>
#include "../include/ParentsImputerByProgeny.h"
#include "../include/common.h"

using namespace std;

ParentsImputerByProgeny::ParentsImputerByProgeny(
							const vector<VCFFamilyRecord *>& rs,
							const vector<vector<int>>& ref_hs_mat,
							const vector<vector<int>>& ref_hs_pat,
							const Map& map_, double w) :
					VCFHMM(map_, w), records(rs),
					ref_haps_mat(ref_hs_mat),
					ref_haps_pat(ref_hs_pat),
					prev_h_table(collect_possible_previous_hidden_states()),
					Cc(calc_Cc(rs)), Cp(calc_Cp(rs)) { }

vector<double> ParentsImputerByProgeny::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> ParentsImputerByProgeny::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

pair<int, int> ParentsImputerByProgeny::compute_parent_gt(int i, int h) const {
	auto parent_gt = [&](int hc, int hp, int prog_a, bool maternal) -> int {
		const auto& ref_haps = maternal ? ref_haps_mat : ref_haps_pat;
		
		if(hc == 0)
			return prog_a | (ref_haps[hp][i] << 1);
		else
			return ref_haps[hp][i] | (prog_a << 1);
	};
	
	const int	hw = h & 1;				// Progeny haplotype inherited from mat
	const int	hc1 = (h >> 1) & 1;		// Maternal haplotype passed to progeny
	const int	hc2 = (h >> 2) & 1;		// Paternal haplotype passed to progeny
	const int	hp1 = (h >> 3) % NH();	// Maternal haplotype not passed on
	const int	hp2 = (h >> 3) / NH();	// Paternal haplotype not passed on
	
	// Allele transmitted from the parents
	const int prog_allele1 = records[i]->get_allele(2, hw);
	const int prog_allele2 = records[i]->get_allele(2, 1 - hw);
	
	return make_pair(
		parent_gt(hc1, hp1, prog_allele1, true),
		parent_gt(hc2, hp2, prog_allele2, false)
	);
}

double ParentsImputerByProgeny::emission_probability(std::size_t i,
															int h) const {
	const auto	p = compute_parent_gt(i, h);
	const int	mat_phased_gt = p.first;
	const int	pat_phased_gt = p.second;
	const int	mat_gt = records[i]->unphased_mat();
	const int	pat_gt = records[i]->unphased_pat();
	return E[mat_phased_gt][mat_gt] + E[pat_phased_gt][pat_gt];
}

double ParentsImputerByProgeny::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cc = Cc[i-1];	// progeny transition probability
	const double	cp = Cp[i-1];	// parental transition probability
	const int	prev_hc1 = (prev_h >> 1) & 1;
	const int	prev_hc2 = (prev_h >> 2) & 1;
	const int	prev_hp1 = (prev_h >> 3) % NH();
	const int	prev_hp2 = (prev_h >> 3) / NH();
	const int	hc1 = (h >> 1) & 1;
	const int	hc2 = (h >> 2) & 1;
	const int	hp1 = (h >> 3) % NH();
	const int	hp2 = (h >> 3) / NH();
	
	return log(prev_hc1 != hc1 ? cc : 1.0 - cc) +
		   log(prev_hc2 != hc2 ? cc : 1.0 - cc) +
		   log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
		   log(prev_hp2 != hp2 ? cp : 1.0 - cp);
}

vector<ParentsImputerByProgeny::DP>
			ParentsImputerByProgeny::initialize_dp() const {
	const size_t	L = num_states();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, h);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ParentsImputerByProgeny::collect_possible_previous_hidden_states() const {
	const int L = num_states();
	vector<vector<int>> prev_h_table(L);
	
	for(int h = 0; h < L; ++h) {
		const int	hw  = h & 1;
		const int	hc1 = (h >> 1) & 1;
		const int	hc2 = (h >> 2) & 1;
		const int	hp1 = (h >> 3) % NH();
		const int	hp2 = (h >> 3) / NH();
		
		// Assume that at most one component can switch simultaneously
		for(int prev_h = hw; prev_h < L; prev_h += 2) {
			const int	prev_hc1 = (prev_h >> 1) & 1;
			const int	prev_hc2 = (prev_h >> 2) & 1;
			const int	prev_hp1 = (prev_h >> 3) % NH();
			const int	prev_hp2 = (prev_h >> 3) / NH();
			
			int counter = static_cast<int>(prev_hc1 != hc1) +
						  static_cast<int>(prev_hc2 != hc2) +
						  static_cast<int>(prev_hp1 != hp1) +
						  static_cast<int>(prev_hp2 != hp2);
			if(counter <= 1)
				prev_h_table[h].push_back(prev_h);
		}
	}
	
	return prev_h_table;
}

void ParentsImputerByProgeny::update_dp(size_t i, vector<DP>& dp) const {
	const int	L = num_states();
	for(int h = 0; h < L; ++h) {	// hidden state
		// Emission log-probability for state h at position i
		const double	E_all = emission_probability(i, h);
		for(int prev_h : prev_h_table[h]) {
			// Transition log-probability from prev_h -> h
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentsImputerByProgeny::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		VCFFamilyRecord	*record = records[i];
		const auto	p = compute_parent_gt(i, hs[i]);
		const int	mat_gt = p.first;
		const int	pat_gt = p.second;
		record->set_geno(0, mat_gt | 4);
		record->set_geno(1, pat_gt | 4);
	}
}

void ParentsImputerByProgeny::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
