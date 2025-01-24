#include <cmath>
#include "../include/VCFNoParentImputed.h"
#include "../include/common.h"

using namespace std;

VCFNoParentImputed::VCFNoParentImputed(
							const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(), VCFFamilyBase(),
				VCFMeasurable(map_), records(rs), ref_haps(ref_hs),
				E{{log(1.0-w*2), log(w/2),     log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(w/2),     log(1.0-w*2), log(w)}},
				Epc{{{log(1.0-w*2),   log(w/2),       log(w/2),       log(w)},
					 {log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)},
					 {log(w/2),       log(1.0-w*2),   log(w/2),       log(w)},
					 {log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)}},
					{{log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)},
					 {log(0.25-w/8),  log(0.5-w*3/4), log(0.25-w/8),  log(w)},
					 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)},
					 {log(0.25+w/4),  log(0.5-w*3/4), log(0.25+w/4),  log(w)}},
					{{log(w/2),       log(1.0-w*2),   log(w/2),       log(w)},
					 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)},
					 {log(w/2),       log(w/2),       log(1.0-w*2),   log(w)},
					 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)}}},
				Cp(calc_Cp(records)) { }

VCFNoParentImputed::~VCFNoParentImputed() {
	Common::delete_all(records);
}

// -> Morgan
double VCFNoParentImputed::dist(const VCFRecord *r1,
								const VCFRecord *r2) const {
	const double	d = (cM(r2->pos()) - cM(r1->pos())) / 100;
	if(d != 0.0)
		return d;
	else	// probably outside map
		return (r2->pos() - r1->pos()) * 1e-6;
}

vector<double> VCFNoParentImputed::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

int VCFNoParentImputed::compute_non_phased_parent_gt(int h, int i) const {
	const size_t	NH = this->ref_haps.size();
	const int	hp1 = h % NH;
	const int	hp2 = h / NH;
	return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
}

double VCFNoParentImputed::emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const size_t	N = num_progenies();
	const auto	*record = records[i];
	const int	mat_phased_gt = compute_non_phased_parent_gt(h, i);
	// emission probability of parent
	const double	Ep = E[mat_phased_gt][mat_gt];
	double	Ec = 0.0;	// emission probability of progenies
	for(size_t j = 0; j < N; ++j) {
		const int	oc = Genotype::gt_to_int(record->get_gt(j+2));
		Ec += Epc[mat_phased_gt][pat_gt][oc];
	}
	return Ep + Ec;
}

vector<VCFNoParentImputed::DP> VCFNoParentImputed::initialize_dp() const {
	const size_t	NH = this->ref_haps.size();
	const size_t	L = NH * NH;
	const size_t	M = this->size();
	vector<DP>	dp(M, DP(L, make_pair(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = records[0];
	const int	mat_gt = Genotype::gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::gt_to_int(record->get_gt(1));
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, h, mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
VCFNoParentImputed::collect_possible_previous_hidden_states() const {
	const size_t	NH = this->ref_haps.size();
	const size_t	L = NH * NH;
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hp1 = h % NH;
		const int	hp2 = h / NH;
		// Neither side will crossover
		for(int h1 = 0; h1 < (int)NH; ++h1) {
			const int	prev_h1 = h1 + hp2 * NH;
			prev_h_table[h].push_back(prev_h1);
		}
		for(int h2 = 0; h2 < (int)NH; ++h2) {
			if(h2 == hp2)	// for duplication
				continue;
			const int	prev_h1 = hp1 + h2 * NH;
			prev_h_table[h].push_back(prev_h1);
		}
	}
	
	return prev_h_table;
}

void VCFNoParentImputed::update_dp(size_t i, vector<DP>& dp,
								const vector<vector<int>>& prev_h_table) const {
	const size_t	NH = this->ref_haps.size();
	const size_t	L = NH * NH;
	// Because the parent is not a direct descendant of the reference,
	// it is easier to cross over than a direct descendant.
	// its multiplication factor
	const VCFFamilyRecord	*record = records[i];
	const double	cp = Cp[i-1];
	
	// observed parent
	const int	mat_gt = Genotype::gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::gt_to_int(record->get_gt(1));
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, mat_gt, pat_gt);
		
		const int	hp1 = h % NH;
		const int	hp2 = h / NH;
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			
			// log of transition probability for parent
			const int	prev_hp1 = prev_h % NH;
			const int	prev_hp2 = prev_h / NH;
			const double	Tp = log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
								 log(prev_hp2 != hp2 ? cp : 1.0 - cp);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

vector<int> VCFNoParentImputed::trace_back(const vector<DP>& dps) const {
	const size_t	M = this->size();
	vector<int>	hs(M, 0);
	
	pair<double, int>	max_pair(-1e300, 0);
	int	max_h = 0;
	const size_t	L = dps[0].size();
	for(int h = 0; h < (int)L; ++h) {
		if(dps.back()[h] > max_pair) {
			max_pair = dps.back()[h];
			max_h = h;
		}
	}
	hs[M-1] = max_h;
	int	prev_h = dps.back()[max_h].second;
	hs[M-2] = prev_h;
	for(int i = (int)M - 2; i > 0; --i) {
		prev_h = dps[i][prev_h].second;
		hs[i-1] = prev_h;
	}
	return hs;
}

void VCFNoParentImputed::set_phased_gt(int gt, VCFFamilyRecord *record) {
	const string	GT = Genotype::int_to_phased_gt(gt);
	record->set_GT(0, GT);
}

void VCFNoParentImputed::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->size(); ++i) {
		const int	h = hs[i];
		VCFFamilyRecord	*record = records[i];
		const int	non_phased_parent_gt = compute_non_phased_parent_gt(h, i);
		set_phased_gt(non_phased_parent_gt, record);
	}
}

void VCFNoParentImputed::impute() {
	const auto	prev_h_table = collect_possible_previous_hidden_states();
	
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->size(); ++i) {
		update_dp(i, dp, prev_h_table);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
