#include <cmath>

using namespace std;

VCFOneParentImputed::VCFOneParentImputed(const std::vector<STRVEC>& header,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>& ref_hs,
							bool is_mat_imputed_, const Map& map_, double w) :
		VCFBase(header), VCFSmallBase(), VCFFamilyBase(), VCFMeasurable(map_),
		records(rs), ref_haps(ref_hs), is_mat_imputed(is_mat_imputed_),
		E{{log(1.0-w*3), log(w),       log(w),       log(w)},
		  {log(w),       log(1.0-w*3), log(w),       log(w)},
		  {log(w),       log(1.0-w*3), log(w),       log(w)},
		  {log(w),       log(w),       log(1.0-w*3), log(w)}} { }

int VCFOneParentImputed::gt_by_haplotypes(int hc1, int hc2,
											int mat_gt, int pat_gt) {
	return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
}

double VCFOneParentImputed::emission_probability(size_t i, int h, int op,
													const vector<int>& ocs,
													phased_parent_gt: int) {
	const size_t	N = num_progenies();
	const int	hp = h >> (N*2);
	const int	non_phased_parent_gt = self.ref_haps[hp1][i] |
										(self.ref_haps[hp2][i] << 1);
	const int	mat_gt = is_mat_imputed ? phased_parent_gt
										: non_phased_parent_gt;
	const int	pat_gt = is_mat_imputed ? non_phased_parent_gt
										: phased_parent_gt;
	// emission probability of parent
	const double	Ep = E[non_phased_parent_gt][op];
	double	Ec = 0.0;	// emission probability of progenies
	for(size_t j = 0; j < N; ++j) {
		const int	hc1 = (h >> (j * 2)) & 1;
		const int	hc2 = (h >> (j * 2 + 1)) & 1;
		const int	gtc = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
		Ec += E[gtc][ocs[j]];
	}
	return Ep + Ec;
}

vector<VCFOneParentImputed::DP> VCFOneParentImputed::initialize_dp(
														size_t L, size_t M) {
	vector<DP>	dp(M, vector(L, pair<double, int>(-1e300, 0)));
	const size_t	phased_col = is_mat_imputed ? 9 : 10;
	non_phased_col = 10 if self.is_mat_imputed else 9
	record = self.records[0]
	phased_parent_gt = phased_gt_to_int(record.v[phased_col])
	op = gt_to_int(record.v[non_phased_col])		# observed parent
	# observed progs
	ocs = [ gt_to_int(gt) for gt in record.v[11:] ]
	for h in range(L):		# hidden state
		E_all = emission_probability(0, h, op, ocs, phased_parent_gt)
		dp[0][h] = (E_all, h)
	return dp


void VCFOneParentImputed::impute() {
	
}
