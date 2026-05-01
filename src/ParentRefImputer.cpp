#include <cmath>

#include "../include/ParentRefImputer.h"

using namespace std;


//////////////////// ParentRefImputer ////////////////////

double ParentRefImputer::progs_emission_probability(std::size_t i,
										int mat_gt, int pat_gt) const {
	const auto	*record = records[i];
	double	Ec = 0.0;
	for(size_t j = 0; j < num_progenies(); ++j) {
		const int	gt = record->get_geno(j+2);
		const int	oc = gt != Genotype::NA ? (gt & 3) : 4;
		Ec += Epc[mat_gt&3][pat_gt&3][oc];
	}
	return Ec;
}

double ParentRefImputer::mat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const int	phased_mat_gt = compute_phased_gt_by_refhaps(h, i);
	const int	mat_gt1 = mat_gt != Genotype::NA ? (mat_gt & 3) : 4;
	// emission probability of parent
	const double	Ep = E[phased_mat_gt][mat_gt1];
	// emission probability of progenies
	double	Ec = progs_emission_probability(i, phased_mat_gt, pat_gt);
	return Ep + Ec;
}

double ParentRefImputer::pat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const int	phased_pat_gt = compute_phased_gt_by_refhaps(h, i);
	const int	pat_gt1 = pat_gt != Genotype::NA ? (pat_gt & 3) : 4;
	// emission probability of parent
	const double	Ep = E[phased_pat_gt][pat_gt1];
	// emission probability of progenies
	double	Ec = progs_emission_probability(i, mat_gt, phased_pat_gt);
	return Ep + Ec;
}

double ParentRefImputer::emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	if(is_mat_ref)
		return pat_emission_probability(i, h, mat_gt, pat_gt);
	else
		return mat_emission_probability(i, h, mat_gt, pat_gt);
}

double ParentRefImputer::transition_probability(size_t i,
													int prev_hp, int hp) const {
	const double	cp = Cp[i-1];
	const int	hp1 = hp % NH();
	const int	hp2 = hp / NH();
	const int	prev_hp1 = prev_hp % NH();
	const int	prev_hp2 = prev_hp / NH();
	return log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
		   log(prev_hp2 != hp2 ? cp : 1.0 - cp);
}

vector<ParentRefImputer::DP> ParentRefImputer::initialize_dp() const {
	const size_t	L = NH() * NH();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = records[0];
	const int	mat_gt = record->mat_gt();
	const int	pat_gt = record->pat_gt();
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, h, mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ParentRefImputer::collect_possible_previous_hidden_states() const {
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

void ParentRefImputer::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = NH() * NH();
	const VCFFamilyRecord	*record = records[i];
	
	// observed parent
	const int	mat_gt = record->mat_gt();
	const int	pat_gt = record->pat_gt();
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E = emission_probability(i, h, mat_gt, pat_gt);
		
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const double	T = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (T + E);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentRefImputer::update_genotypes(const vector<int>& hs) {
	const size_t	j = phased_index();
	for(size_t i = 0; i < this->M(); ++i) {
		const int	phased_gt = compute_phased_gt_by_refhaps(hs[i], i);
		records[i]->set_geno(j, phased_gt | 4);
	}
}

void ParentRefImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}

array<array<array<double, 5>, 4>, 4> ParentRefImputer::calc_Epc(double w) {
	// hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
	// observed 0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3, N/A: 4
	// Assume a 0.1 probability of confusing 0|1 and 1|0.
	array<array<array<double,5>,4>,4> E1{};
	
	for(int i = 0; i < 16; ++i) {
		const int	mat = i >> 2;
		const int	pat = i & 3;
		array<double, 5>	orig_prob{};
		for(int j = 0; j < 4; ++j) {
			const int	j_mat = j >> 1;
			const int	j_pat = j & 1;
			const int	a1 = (mat >> j_mat) & 1;
			const int	a2 = (pat >> j_pat) & 1;
			const int	gt = a1 | (a2 << 1);
			if(gt == 1) {
				orig_prob[1] += 0.25 * 0.9;
				orig_prob[2] += 0.25 * 0.1;
			}
			else if(gt == 2) {
				orig_prob[1] += 0.25 * 0.1;
				orig_prob[2] += 0.25 * 0.9;
			}
			else {
				orig_prob[gt] += 0.25;
			}
		}
		
		// put errors
		for(int gt = 0; gt < 4; ++gt) {
			for(int error_gt = 0; error_gt < 5; ++error_gt) {
				if(error_gt == 4)
					E1[mat][pat][error_gt] = w;
				else if(error_gt == gt)
					E1[mat][pat][error_gt] += orig_prob[gt] * (1.0 - w*2);
				else
					E1[mat][pat][error_gt] += orig_prob[gt] * (w/3);
			}
		}
	}
	
	array<array<array<double,5>,4>,4> E{};
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			for(int k = 0; k < 5; ++k)
				E[i][j][k] = log(E1[i][j][k]);
		}
	}
	return E;
}

vector<double> ParentRefImputer::calc_Cp() const {
	const double	K = 5.0;
	const size_t	M = records.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(records[i], records[i+1]) * K);
	}
	return Cp;
}
