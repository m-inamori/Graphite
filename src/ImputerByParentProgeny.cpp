#include <cmath>
#include "../include/ImputerByParentProgeny.h"
#include "../include/common.h"

using namespace std;

ImputerByParentProgeny::ImputerByParentProgeny(
									const vector<VCFFamilyRecord *>& rs,
									const vector<vector<int>>& ref_hs,
									bool should_impute_mat_,
									const Map& map_, double w) :
					VCFHMM(map_, w), records(rs), ref_haps(ref_hs),
					prev_h_table(collect_possible_previous_hidden_states()),
					should_impute_mat(should_impute_mat_),
					Cc(calc_Cc(rs)), Cp(calc_Cp(rs)) { }

vector<double> ImputerByParentProgeny::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> ImputerByParentProgeny::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

void ImputerByParentProgeny::decode(int h, int& hw, int& hp1, int& hp21,
											int& hp22, vector<int>& hc) const {
	hw = h & 1;
	hp1 = (h >> 1) & 1;
	hp21 = (h >> 2) & 1;
	const size_t	np = num_should_impute_progenies();
	for(size_t i = 3; i < np*2+3; ++i) {
		hc.push_back((h >> i) & 1);
	}
	hp22 = h >> (np*2 + 3);
}

int ImputerByParentProgeny::parent_gt_by_haplotypes(size_t i, int hw,
													int hp1, int hp2,
													int prog_gt) const {
	if(hp1 == 0)
		return ((prog_gt >> (1 - hw)) & 1) | (ref_haps[hp2][i] << 1);
	else
		return (ref_haps[hp2][i] << 1) | ((prog_gt >> (1 - hw)) & 1);
}

int ImputerByParentProgeny::progeny_gt_by_haplotypes(size_t j,
											const std::vector<int>& hc,
											int mat_gt, int pat_gt) const {
	return ((mat_gt >> hc[j * 2]) & 1) | ((pat_gt >> hc[j * 2 + 1]) & 1);
}

vector<int> ImputerByParentProgeny::compute_gts(size_t i, int h) const {
	std::vector<int> gts(num_should_impute_progenies() + 1, 0);
	int	hw, hp1, hp21, hp22;
	vector<int>	hc;
	decode(h, hw, hp1, hp21, hp22, hc);
	
	const auto	*record = records[i];
	const int	prog_gt = record->get_geno(2);
	const int parent_gt = parent_gt_by_haplotypes(i, hw, hp21, hp22, prog_gt);
	gts[0] = parent_gt;
	const int	mat_gt = should_impute_mat ? parent_gt : record->get_geno(0);
	const int	pat_gt = should_impute_mat ? record->get_geno(1) : parent_gt;
	for(size_t j = 0; j < num_should_impute_progenies(); ++j) {
		gts[j+1] = progeny_gt_by_haplotypes(j, hc, mat_gt, pat_gt);
	}
	return gts;
}

double ImputerByParentProgeny::emission_probability(size_t i, int h) const {
	const auto	*record = records[i];
	const vector<int>	phased_gts = compute_gts(i, h);
	int	hw, hp1, hp21, hp22;
	vector<int>	hc;
	decode(h, hw, hp1, hp21, hp22, hc);
	
	const size_t	index = should_imputed_index();
	
	// Check whether the allele transmitted from the imputed parent
	// to the imputed progeny is consistent.
	const int	parent_a = (record->get_geno(1-index) >> hp1) & 1;
	const int	progeny_a = (record->get_geno(2) >> hw) & 1;
	double	E = (parent_a == progeny_a) ? std::log(0.99) : std::log(0.01);
	const int parent_gt = record->unphased(index);
	E += this->E[phased_gts[0]][parent_gt];
	for(size_t j = 0; j < num_should_impute_progenies(); ++j) {
		const int	gt = record->unphased(j + 3);
		E += this->E[phased_gts[j+1]][gt];
	}
	return E;
}

double ImputerByParentProgeny::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cc = Cc[i - 1];	// Progeny transition probability
	const double	cp = Cp[i - 1];	// Parent transition probability
	
	int	prev_hw, prev_hp1, prev_hp21, prev_hp22;
	vector<int>	prev_hc;
	decode(prev_h, prev_hw, prev_hp1, prev_hp21, prev_hp22, prev_hc);
	
	int	hw, hp1, hp21, hp22;
	vector<int>	hc;
	decode(h, hw, hp1, hp21, hp22, hc);
	
	// Transition probability:
	// use c when the state changes, otherwise use (1 - c).
	double	T = log(prev_hp1 != hp1 ? cc : 1.0 - cc) +
				log(prev_hp21 != hp21 ? cc : 1.0 - cc) +
				log(prev_hp22 != hp22 ? cp : 1.0 - cp);
	for(size_t j = 0; j < hc.size(); ++j) {
		T += log(prev_hc[j] != hc[j] ? cc : 1.0 - cc);
	}
	return T;
}

vector<ImputerByParentProgeny::DP>
					ImputerByParentProgeny::initialize_dp() const {
	const size_t	L = num_states();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	for(int h = 0; h < static_cast<int>(L); ++h) {
		const double	E_all = emission_probability(0, h);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ImputerByParentProgeny::collect_possible_previous_hidden_states() const {
	const size_t	L = num_states();
	vector<vector<int>>	prev_h_table(L);
	for(int h = 0; h < static_cast<int>(L); ++h) {
		int	hw, hp1, hp21, hp22;
		vector<int>	hc;
		decode(h, hw, hp1, hp21, hp22, hc);
		
		// Assume that multiple crossover events cannot occur simultaneously.
		for(int prev_h = hw; prev_h < static_cast<int>(L); prev_h += 2) {
			int	prev_hw, prev_hp1, prev_hp21, prev_hp22;
			vector<int>	prev_hc;
			decode(prev_h, prev_hw, prev_hp1, prev_hp21, prev_hp22, prev_hc);
			
			int	counter = static_cast<int>(prev_hp1 != hp1) +
						  static_cast<int>(prev_hp21 != hp21) +
						  static_cast<int>(prev_hp22 != hp22);
			for(size_t j = 0; j < hc.size(); ++j) {
				counter += static_cast<int>(prev_hc[j] != hc[j]);
			}
			
			if(counter <= 1)
				prev_h_table[h].push_back(prev_h);
		}
	}
	
	return prev_h_table;
}

void ImputerByParentProgeny::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = num_states();
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

void ImputerByParentProgeny::update_genotypes(const vector<int>& hs) {
	const size_t	index = should_imputed_index();
	for(size_t i = 0; i < this->M(); ++i) {
		auto	*record = records[i];
		const vector<int>	gts = compute_gts(i, hs[i]);
		record->set_geno(index, gts[0] | 4);
		for(size_t j = 3; j < num_should_impute_progenies() + 3; ++j) {
			record->set_geno(j, gts[j-2] | 4);
		}
	}
}

void ImputerByParentProgeny::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
