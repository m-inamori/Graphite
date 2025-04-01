#include <cmath>
#include "../include/ParentImputerByProgeny.h"
#include "../include/common.h"

using namespace std;

ParentImputerByProgeny::ParentImputerByProgeny(const vector<VCFRecord *>& rs,
							 const vector<vector<int>>& ref_hs,
							 bool is_mat_known_, const Map& map_, double w) :
					VCFHMM(rs, map_, w), ref_records(rs), ref_haps(ref_hs),
					prev_h_table(collect_possible_previous_hidden_states()),
					is_mat_known(is_mat_known_),
					Cc(calc_Cc(rs)), Cp(calc_Cp(rs)) { }

vector<double> ParentImputerByProgeny::calc_Cc(
							const vector<VCFRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> ParentImputerByProgeny::calc_Cp(
							const vector<VCFRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

int ParentImputerByProgeny::compute_parent_gt(size_t i, int h) const {
	const int	hw = h & 1;			// Which of the two haplotypes was passed
									// from the parent to the progeny
	const int	hc = (h >> 1) & 1;	// Which parent passed it to the progeny
	const int	hp = h >> 2;		// Those who have not passed to the progeny
	const string&	prog_gt = records[i]->get_gt(1);
	const int	prog_allele = static_cast<int>(prog_gt.c_str()[hw*2] - '0');
	if(hc == 0)
		return prog_allele | (ref_haps[hp][i] << 1);
	else
		return ref_haps[hp][i] | (prog_allele << 1);
}

double ParentImputerByProgeny::emission_probability(size_t i, int h) const {
	const int	parent_phased_gt = compute_parent_gt(i, h);
	const int	parent_gt = Genotype::gt_to_int(records[i]->get_GT(0));
	return E[parent_phased_gt][parent_gt];
}

double ParentImputerByProgeny::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cc = Cc[i-1];
	const double	cp = Cp[i-1];
	const int	hc = (h >> 1) & 1;
	const int	hp = h >> 2;
	const int	prev_hc = (prev_h >> 1) & 1;
	const int	prev_hp = prev_h >> 2;
	return log(prev_hc != hc ? cc : 1.0 - cc) +
		   log(prev_hp != hp ? cp : 1.0 - cp);
}

vector<ParentImputerByProgeny::DP>
					ParentImputerByProgeny::initialize_dp() const {
	const size_t	L = NH() << 2;
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(0, h);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ParentImputerByProgeny::collect_possible_previous_hidden_states() const {
	const size_t	L = NH() << 2;
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hw = h & 1;		// invariant
		const int	hc = (h >> 1) & 1;
		const int	hp = h >> 2;
		// Neither side will crossover
		for(int prev_h = hw; prev_h < (int)L; prev_h += 2) {
			const int	prev_hc = (prev_h >> 1) & 1;
			const int	prev_hp = prev_h >> 2;
			if(prev_hp == hp || prev_hc == hc)
				prev_h_table[h].push_back(prev_h);
		}
	}
	
	return prev_h_table;
}

void ParentImputerByProgeny::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = NH() << 2;
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h);
		for(auto p = prev_h_table[h].begin(); p != prev_h_table[h].end(); ++p) {
			const int	prev_h = *p;
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentImputerByProgeny::update_genotypes(const vector<int>& hs) {
	// Genotypes may be swapped depending on
	// which of the progeny passed the parental haplotype
	const bool	is_swapped = is_mat_known ^ ((hs[0] & 1) == 0);
	for(size_t i = 0; i < this->M(); ++i) {
		VCFRecord	*record = ref_records[i];
		const int	parent_gt = compute_parent_gt(i, hs[i]);
		record->set_GT(0, Genotype::int_to_phased_gt(parent_gt));
		if(is_swapped) {
			string&	gt = record->get_mut_gt(0);
			std::swap(gt[0], gt[2]);
		}
	}
}

void ParentImputerByProgeny::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
