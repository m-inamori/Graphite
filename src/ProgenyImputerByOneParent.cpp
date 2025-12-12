#include <cmath>
#include "../include/ProgenyImputerByOneParent.h"
#include "../include/common.h"

using namespace std;

ProgenyImputerByOneParent::ProgenyImputerByOneParent(
											const vector<VCFFamilyRecord *>& rs,
											const vector<vector<int>>& ref_hs,
											bool is_mat,
											const Map& map_, double w) :
										VCFHMM(rs, map_, w), ref_records(rs),
										ref_haps(ref_hs),
										is_mat_imputed(is_mat),
										Cc(calc_Cc(rs)), Cp(calc_Cp(rs)) { }

vector<double> ProgenyImputerByOneParent::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> ProgenyImputerByOneParent::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

double ProgenyImputerByOneParent::emission_probability(size_t i, size_t j,
														int h,
														int parent_gt) const {
	const int	oc = ref_records[i]->unphased(j+2);
	const auto	alleles = parent_alleles(h, i);
	const int	mat_a = alleles.first;
	const int	pat_a = alleles.second;
	const int	phased_gt = Genotype::from_alleles(mat_a, pat_a) & 3;
	return E[phased_gt][oc];
}

vector<ProgenyImputerByOneParent::DP>
ProgenyImputerByOneParent::initialize_dp(size_t j) const {
	const size_t	L = NH() * 2;
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = ref_records[0];
	const int	parent_gt = is_mat_imputed ? record->mat_gt() :
												record->pat_gt();
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, j, h, parent_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void ProgenyImputerByOneParent::update_dp(size_t i, size_t j,
												vector<DP>& dp) const {
	const size_t	L = NH() * 2;
	const VCFFamilyRecord	*record = ref_records[i];
	const int	parent_gt = is_mat_imputed ? record->mat_gt() :
												record->pat_gt();
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E = emission_probability(i, j, h, parent_gt);
		for(int prev_h = 0; prev_h < (int)L; ++prev_h) {
			const double	T = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (T + E);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ProgenyImputerByOneParent::update_genotypes(size_t j,
													const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		const auto	alleles = parent_alleles(hs[i], i);
		const int	mat_a = alleles.first;
		const int	pat_a = alleles.second;
		const int	phased_gt = Genotype::from_alleles(mat_a, pat_a);
		ref_records[i]->set_geno(j+2, phased_gt);
	}
}

void ProgenyImputerByOneParent::impute(size_t j) {
	// DP
	vector<DP>	dp = initialize_dp(j);
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, j, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(j, hs);
}
