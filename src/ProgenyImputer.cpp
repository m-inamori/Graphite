#include <cmath>
#include "../include/ProgenyImputer.h"
#include "../include/common.h"

using namespace std;

ProgenyImputer::ProgenyImputer(const vector<VCFFamilyRecord *>& rs,
												const Map& map_, double w) :
				VCFHMM(rs, map_, w), ref_records(rs), Cc(calc_Cc(rs)) { }

vector<double> ProgenyImputer::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

double ProgenyImputer::emission_probability(size_t i, size_t j, int h,
												int mat_gt, int pat_gt) const {
	const auto	*record = ref_records[i];
	const int	oc = Genotype::gt_to_int(record->get_gt(j+2));
	const int	phased_gt = gt_by_haplotypes(h, mat_gt, pat_gt);
	return E[phased_gt][oc];
}

vector<ProgenyImputer::DP> ProgenyImputer::initialize_dp(size_t j) const {
	vector<DP>	dp(M(), DP(4, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = ref_records[0];
	const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
	for(int h = 0; h < 4; ++h) {	// hidden state
		const double	E_all = emission_probability(0, j, h, mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void ProgenyImputer::update_dp(size_t i, size_t j, vector<DP>& dp) const {
	const VCFFamilyRecord	*record = ref_records[i];
	
	// observed parent
	const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
	
	for(int h = 0; h < 4; ++h) {
		const double	E_all = emission_probability(i, j, h, mat_gt, pat_gt);
		for(int prev_h = 0; prev_h < 4; ++prev_h) {
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ProgenyImputer::update_genotypes(size_t j, const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		const VCFFamilyRecord	*record = ref_records[i];
		const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
		const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
		const int	phased_gt = gt_by_haplotypes(hs[i], mat_gt, pat_gt);
		ref_records[i]->set_GT(j+2, Genotype::int_to_phased_gt(phased_gt));
	}
}

void ProgenyImputer::impute(size_t j) {
	// DP
	vector<DP>	dp = initialize_dp(j);
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, j, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(j, hs);
}
