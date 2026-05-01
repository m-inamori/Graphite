#include <cmath>
#include "../include/ProgenySelfRefImputer.h"
#include "../include/common.h"

using namespace std;


//////////////////// ProgenySelfRefImputer ////////////////////

ProgenySelfRefImputer::ProgenySelfRefImputer(const vector<GenoRecord *>& rs,
												const Map& map_, double w) :
									VCFHMMSelfRef(map_, w),
									records(rs),
									Cc(calc_Cc(rs)) { }

vector<double> ProgenySelfRefImputer::calc_Cc(
							const vector<GenoRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

double ProgenySelfRefImputer::emission_probability(size_t i, size_t j, int h,
														int parent_gt) const {
	const int	gt = records[i]->get_geno(j+1);
	const int	oc = gt == Genotype::NA ? 4 : gt & 3;
	const int	phased_gt = gt_by_haplotypes(h, parent_gt);
	return E[phased_gt][oc];
}

vector<VCFHMMSelfRef::DP> ProgenySelfRefImputer::initialize_dp(size_t j) const {
	vector<DP>	dp(M(), DP(4, pair<double, int>(MIN_PROB, 0)));
	const GenoRecord	*record = records[0];
	const int	parent_gt = record->get_geno(0);
	for(int h = 0; h < 4; ++h) {	// hidden state
		const double	E = emission_probability(0, j, h, parent_gt);
		dp[0][h] = make_pair(E, h);
	}
	return dp;
}

void ProgenySelfRefImputer::update_dp(size_t i, size_t j,
											vector<DP>& dp) const {
	const GenoRecord	*record = records[i];
	const int	parent_gt = record->get_geno(0) & 3;
	
	for(int h = 0; h < 4; ++h) {
		const double	E = emission_probability(i, j, h, parent_gt);
		
		for(int prev_h = 0; prev_h < 4; ++prev_h) {
			const double	T = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (T + E);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ProgenySelfRefImputer::update_genotypes(size_t j, const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		const GenoRecord	*record = records[i];
		const int	parent_gt = record->get_geno(0);
		const int	gtc_int = gt_by_haplotypes(hs[i], parent_gt);
		records[i]->set_geno(j+1, gtc_int | 4);
	}
}

void ProgenySelfRefImputer::impute(size_t j) {
	// DP
	vector<DP>	dp = initialize_dp(j);
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, j, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(j, hs);
}
