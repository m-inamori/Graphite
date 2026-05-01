#include <cmath>
#include "../include/ParentNoRefImputer.h"
#include "../include/common.h"

using namespace std;

vector<double> ParentNoRefImputer::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

double ParentNoRefImputer::progs_emission_probability(size_t i, int mat_gt,
															int pat_gt) const {
	double	Ec = 0.0;
	for(size_t j = 0; j < num_progenies(); ++j) {
		const int	gt = records[i]->get_geno(j+2);
		const int	oc = gt != Genotype::NA ? gt & 3 : 4;
		Ec += Epc[mat_gt&3][pat_gt&3][oc];
	}
	return Ec;
}

double ParentNoRefImputer::mat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const{
	const int	phased_mat_gt = compute_phased_gt_by_refhaps(h, i);
	const int	mat_gt1 = mat_gt != Genotype::NA ? mat_gt & 3 : 4;
	const double	Ep = E[phased_mat_gt][mat_gt1];		// mat emission
	// probability of progenies' emission
	const double	Ec = progs_emission_probability(i, phased_mat_gt, pat_gt);
	return Ep + Ec;
}

double ParentNoRefImputer::pat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const int	phased_pat_gt = compute_phased_gt_by_refhaps(h, i);
	const int	pat_gt1 = pat_gt != Genotype::NA ? pat_gt & 3 : 4;
	const double	Ep = E[phased_pat_gt][pat_gt1];		// pat emission
	// probability of progenies' emission
	const double	Ec = progs_emission_probability(i, mat_gt, phased_pat_gt);
	return Ep + Ec;
}

double ParentNoRefImputer::emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	if(is_mat)
		return mat_emission_probability(i, h, mat_gt, pat_gt);
	else
		return pat_emission_probability(i, h, mat_gt, pat_gt);
}

double ParentNoRefImputer::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cp = Cp[i-1];
	const int	h1 = h % NH();
	const int	h2 = h / NH();
	const int	prev_h1 = prev_h % NH();
	const int	prev_h2 = prev_h / NH();
	return log(prev_h1 != h1 ? cp : 1.0 - cp) +
		   log(prev_h2 != h2 ? cp : 1.0 - cp);
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ParentNoRefImputer::collect_possible_previous_hidden_states() const {
	vector<vector<int>>	prev_h_table(L(), vector<int>());
	for(int h = 0; h < (int)L(); ++h) {		// hidden state
		const int	hp1 = h % NH();
		const int	hp2 = h / NH();
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

vector<ParentNoRefImputer::DP>
					ParentNoRefImputer::initialize_dp() const {
	vector<DP>	dp(M(), DP(L(), pair<double, int>(MIN_PROB, 0)));
	const int	mat_gt = records[0]->mat_gt();
	const int	pat_gt = records[0]->pat_gt();
	for(int h = 0; h < (int)L(); ++h) {
		const double	E_all = emission_probability(0, h, mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void ParentNoRefImputer::update_dp(size_t i, vector<DP>& dp) const {
	const int	mat_gt = records[i]->mat_gt();
	const int	pat_gt = records[i]->pat_gt();
	for(int h = 0; h < (int)L(); ++h) {
		const double	E_all = emission_probability(i, h, mat_gt, pat_gt);
		for(auto p = prev_h_table[h].begin(); p != prev_h_table[h].end(); ++p) {
			const int	prev_h = *p;
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentNoRefImputer::update_genotypes(const vector<int>& hs) {
	// Genotypes may be swapped depending on
	// which of the progeny passed the parental haplotype
	const size_t	j = is_mat ? 0 : 1;
	for(size_t i = 0; i < M(); ++i) {
		const int	parent_gt = compute_phased_gt_by_refhaps(hs[i], i);
		records[i]->set_geno(j, parent_gt | 4);
	}
}

void ParentNoRefImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}

ParentNoRefImputer::EpcArray ParentNoRefImputer::calc_Epc(double w) {
	auto divide_into_alleles = [](int gt) -> std::array<int,2> {
		if(gt == 4) {	// NA
			return {0, 1};
		}
		else {
			return {gt & 1, gt >> 1};
		}
	};
	
	// hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
	// observed 0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3, N/A: 4
	// Assume that 0|1 <=> 1|0 is incorrect with a probability of 0.1.
	EpcArray	E1{};	// all elements initialize to 0.0
	for(int index = 0; index < 25; ++index) {
		const int	mat = index / 5;
		const int	pat = index % 5;
		array<double, 5>	orig_prob{};
		const auto	as1 = divide_into_alleles(mat);
		const auto	as2 = divide_into_alleles(pat);
		for(int i = 0; i < 4; ++i) {
			const int	a1 = as1[i>>1];
			const int	a2 = as2[i&1];
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
	
	EpcArray E{};
	for(int i = 0; i < 5; ++i) {
		for(int j = 0; j < 5; ++j) {
			for(int k = 0; k < 5; ++k)
				E[i][j][k] = log(E1[i][j][k]);
		}
	}
	return E;
}
