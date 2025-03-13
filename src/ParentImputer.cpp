#include <cmath>
#include "../include/ParentImputer.h"
#include "../include/common.h"

using namespace std;

ParentImputer::ParentImputer(const vector<VCFFamilyRecord *>& rs,
							 bool is_mat_imputed_,
							 const vector<vector<int>>& ref_hs,
							 const Map& map_, double w) :
			VCFHMM(rs, map_, w), ref_records(rs), ref_haps(ref_hs),
			prev_h_table(collect_possible_previous_hidden_states()),
			is_mat_imputed(is_mat_imputed_), Cp(calc_Cp(rs)),
			Epc{{{log(1.0-w*2),   log(w/2),       log(w/2),       log(w)},
				 {log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)},
				 {log(w/2),       log(1.0-w*2),   log(w/2),       log(w)},
				 {log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)}},
				{{log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)},
				 {log(0.25-w/8),  log(0.5-w*3/4), log(0.25-w/8),  log(w)},
				 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)},
				 {log(0.25-w/8),  log(0.5-w*3/4), log(0.25-w/8),  log(w)}},
				{{log(w/2),       log(1.0-w*2),   log(w/2),       log(w)},
				 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)},
				 {log(w/2),       log(w/2),       log(1.0-w*2),   log(w)},
				 {log(w/2),       log(0.5-w*3/4), log(0.5-w*3/4), log(w)}},
				{{log(0.5-w*3/4), log(0.5-w*3/4), log(w/2),       log(w)},
				 {log(0.25-w/8),  log(0.5-w*3/4), log(0.25-w/8),  log(w)},
				  {log(w/2),      log(0.5-w*3/4), log(0.5-w*3/4), log(w)},
				  {log(0.25),     log(0.25),      log(0.25),      log(0.25)}}}
			{ }

vector<double> ParentImputer::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

double ParentImputer::progs_emission_probability(std::size_t i,
										int mat_gt, int pat_gt) const {
	const auto	*record = ref_records[i];
	double	Ec = 0.0;
	for(size_t j = 0; j < num_progenies(); ++j) {
		const int	oc = Genotype::gt_to_int(record->get_gt(j+2));
		Ec += Epc[mat_gt][pat_gt][oc];
	}
	return Ec;
}

double ParentImputer::mat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const int	phased_mat_gt = compute_phased_gt_by_refhaps(h, i);
	// emission probability of parent
	const double	Ep = E[phased_mat_gt][mat_gt];
	const int	non_phased_mat_gt = (phased_mat_gt & 1) + (phased_mat_gt >> 1);
	// emission probability of progenies
	double	Ec = progs_emission_probability(i, non_phased_mat_gt, pat_gt);
	return Ep + Ec;
}

double ParentImputer::pat_emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	const int	phased_pat_gt = compute_phased_gt_by_refhaps(h, i);
	// emission probability of parent
	const double	Ep = E[phased_pat_gt][pat_gt];
	const int	non_phased_pat_gt = (phased_pat_gt & 1) + (phased_pat_gt >> 1);
	// emission probability of progenies
	double	Ec = progs_emission_probability(i, mat_gt, non_phased_pat_gt);
	return Ep + Ec;
}

double ParentImputer::emission_probability(size_t i, int h,
												int mat_gt, int pat_gt) const {
	if(is_mat_imputed)
		return mat_emission_probability(i, h, mat_gt, pat_gt);
	else
		return pat_emission_probability(i, h, mat_gt, pat_gt);
}

double ParentImputer::parent_transition_probability(size_t i,
													int prev_hp, int hp) const {
	const double	cp = Cp[i-1];
	const int	hp1 = hp % NH();
	const int	hp2 = hp / NH();
	const int	prev_hp1 = prev_hp % NH();
	const int	prev_hp2 = prev_hp / NH();
	return log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
		   log(prev_hp2 != hp2 ? cp : 1.0 - cp);
}

vector<ParentImputer::DP> ParentImputer::initialize_dp() const {
	const size_t	L = NH() * NH();
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = ref_records[0];
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
ParentImputer::collect_possible_previous_hidden_states() const {
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

void ParentImputer::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	L = NH() * NH();
	const VCFFamilyRecord	*record = ref_records[i];
	
	// observed parent
	const int	mat_gt = Genotype::gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::gt_to_int(record->get_gt(1));
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, mat_gt, pat_gt);
		
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const double	Tp = parent_transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentImputer::update_genotypes(const vector<int>& hs) {
	const size_t	c = phased_col();
	for(size_t i = 0; i < this->M(); ++i) {
		const int	phased_gt = compute_phased_gt_by_refhaps(hs[i], i);
		ref_records[i]->set_GT(c-9, Genotype::int_to_phased_gt(phased_gt));
	}
}

void ParentImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
