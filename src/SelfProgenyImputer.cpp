#include <cmath>
#include "../include/SelfProgenyImputer.h"
#include "../include/common.h"

using namespace std;

SelfProgenyImputer::SelfProgenyImputer(const vector<GenoRecord *>& rs,
									 size_t iprog, const Map& map_, double w) :
					VCFHMM(rs, map_, w), records(rs), ic(iprog), Cc(calc_Cc(rs))
					{ }

vector<double> SelfProgenyImputer::calc_Cc(
							const vector<GenoRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

double SelfProgenyImputer::emission_probability(int h, int parent_gt,
																int oc) const {
	const int	gt_prog = progeny_genotype(h, parent_gt);
	return E[gt_prog][oc];
}

double SelfProgenyImputer::transition_probability(size_t i,
													int prev_h, int h) const {
	const double	cc = Cc[i-1];
	const int	d = h ^ prev_h;
	return (log((d & 1) == 1 ? cc : 1.0 - cc) +
			log((d >> 1) == 1 ? cc : 1.0 - cc));
}

vector<SelfProgenyImputer::DP> SelfProgenyImputer::initialize_dp() const {
	vector<DP>	dp(M(), DP(4, pair<double, int>(MIN_PROB, 0)));
	const GenoRecord	*record = records[0];
	const int	parent_gt = record->get_geno()[0];
	const int	oc = record->get_geno()[ic+1];
	for(int h = 0; h < 4; ++h) {
		const double	E_all = emission_probability(h, parent_gt, oc);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void SelfProgenyImputer::update_dp(size_t i, vector<DP>& dp) const {
	const GenoRecord	*record = records[i];
	const int	parent_gt = record->get_geno()[0];
	const int	oc = record->get_geno()[ic+1];
	for(int h = 0; h < 4; ++h) {
		const double	E_all = emission_probability(h, parent_gt, oc);
		
		for(int prev_h = 0; prev_h < 4; ++prev_h) {
			const double	Tp = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void SelfProgenyImputer::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->M(); ++i) {
		GenoRecord	*record = records[i];
		const int	parent_gt = record->get_geno()[0];
		const int	prog_gt = progeny_genotype(hs[i], parent_gt);
		record->set_geno(ic + 1, prog_gt | 4);
	}
}

void SelfProgenyImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
