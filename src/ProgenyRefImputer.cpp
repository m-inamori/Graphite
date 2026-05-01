#include <cmath>
#include "../include/ProgenyRefImputer.h"
#include "../include/common.h"

using namespace std;

vector<double> ProgenyRefImputer::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

int ProgenyRefImputer::gt_by_haplotypes(int hc, int mat_gt, int pat_gt) const {
	const int	hc1 = hc & 1;
	const int	hc2 = hc >> 1;
	return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
}

double ProgenyRefImputer::emission_probability(size_t i, size_t j, int h,
												int mat_gt, int pat_gt) const {
	const auto	*record = records[i];
	const int	gt = record->get_geno(j+2);
	const int	oc = gt == Genotype::NA ? 4 : (gt & 3);
	const int	phased_gt = gt_by_haplotypes(h, mat_gt, pat_gt);
	return E[phased_gt][oc];
}

double ProgenyRefImputer::progeny_transition_probability(size_t i,
													int prev_hc, int hc) const {
	const double	cc = Cc[i-1];
	const int	hc1 = hc & 1;
	const int	hc2 = hc >> 1;
	const int	prev_hc1 = prev_hc & 1;
	const int	prev_hc2 = prev_hc >> 1;
	return log(prev_hc1 != hc1 ? cc : 1.0 - cc) +
		   log(prev_hc2 != hc2 ? cc : 1.0 - cc);
}

vector<ProgenyRefImputer::DP> ProgenyRefImputer::initialize_dp(size_t j) const {
	vector<DP>	dp(M(), DP(4, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = records[0];
	const int	mat_gt = record->mat_gt();
	const int	pat_gt = record->pat_gt();
	for(int h = 0; h < 4; ++h) {	// hidden state
		const double	E_all = emission_probability(0, j, h, mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void ProgenyRefImputer::update_dp(size_t i, size_t j, vector<DP>& dp) const {
	const VCFFamilyRecord	*record = records[i];
	const int	mat_gt = record->mat_gt() & 3;
	const int	pat_gt = record->pat_gt() & 3;
	
	for(int hc = 0; hc < 4; ++hc) {
		const double	E_all = emission_probability(i, j, hc, mat_gt, pat_gt);
		
		for(int prev_hc = 0; prev_hc < 4; ++prev_hc) {
			const double	T = progeny_transition_probability(i, prev_hc, hc);
			const double	prob = dp[i-1][prev_hc].first + (T + E_all);
			dp[i][hc] = max(dp[i][hc], make_pair(prob, prev_hc));
		}
	}
}

void ProgenyRefImputer::update_genotypes(size_t j, const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		VCFFamilyRecord	*record = records[i];
		const int	mat_gt = record->mat_gt();
		const int	pat_gt = record->pat_gt();
		const int	gtc_int = gt_by_haplotypes(hs[i], mat_gt, pat_gt);
		record->set_geno(j+2, gtc_int | 4);
	}
}

void ProgenyRefImputer::impute(size_t j) {
	// DP
	vector<DP>	dp = initialize_dp(j);
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, j, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(j, hs);
}
