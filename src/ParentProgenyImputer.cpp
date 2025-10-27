#include <cmath>
#include "../include/ParentProgenyImputer.h"
#include "../include/common.h"

using namespace std;

ParentProgenyImputer::ParentProgenyImputer(const vector<VCFFamilyRecord *>& rs,
							 bool is_mat_imputed_,
							 const vector<vector<int>>& ref_hs,
							 const Map& map_, double w) :
					VCFHMM(rs, map_, w), ref_records(rs), ref_haps(ref_hs),
					prev_h_table(collect_possible_previous_hidden_states()),
					is_mat_imputed(is_mat_imputed_),
					Cc(calc_Cc(rs)), Cp(calc_Cp(rs)) { }

vector<double> ParentProgenyImputer::calc_Cc(
							const vector<VCFFamilyRecord *>& rs) const {
	const size_t	M = rs.size();
	vector<double>	Cc(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cc[i] = Map::Kosambi(dist(rs[i], rs[i+1]));
	}
	return Cc;
}

vector<double> ParentProgenyImputer::calc_Cp(
							const vector<VCFFamilyRecord *>& rs) const {
	const double	K = 5.0;
	const size_t	M = rs.size();
	vector<double>	Cp(M-1);
	for(size_t i = 0; i < M - 1; ++i) {
		Cp[i] = Map::Kosambi(dist(rs[i], rs[i+1]) * K);
	}
	return Cp;
}

double ParentProgenyImputer::emission_probability(std::size_t i, int h, int op,
												const std::vector<int>& ocs,
												int phased_parent_gt) const {
	const size_t	N = num_progenies();
	const int	hp = h >> (N*2);
	const int	hp1 = hp % NH();
	const int	hp2 = hp / NH();
	const int	non_phased_parent_gt = ref_haps[hp1][i] |
										ref_haps[hp2][i] << 1;
	const int	mat_gt = is_mat_imputed ? phased_parent_gt :
												non_phased_parent_gt;
	const int	pat_gt = is_mat_imputed ? non_phased_parent_gt :
												phased_parent_gt;
	const double	Ep = E[non_phased_parent_gt][op];
	double	Ec = 0.0;
	for(size_t j = 0; j < N; ++j) {
		const int	hc = (h >> (j * 2)) & 3;
		const int	hc1 = hc & 1;
		const int	hc2 = hc >> 1;
		const int	gtc = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt);
		Ec += E[gtc][ocs[j]];
	}
	return Ep + Ec;
}

double ParentProgenyImputer::progeny_transition_probability(size_t i,
													int prev_hc, int hc) const {
	const int	t = prev_hc ^ hc;
	const double	cc = Cc[i-1];
	double	Tc = 0.0;
	for(size_t j = 0; j < num_progenies() * 2; ++j) {
		const int	tc = (t >> j) & 1;
		Tc += log(tc != 0 ? cc : 1.0 - cc);
	}
	return Tc;
}

double ParentProgenyImputer::parent_transition_probability(size_t i,
													int prev_hp, int hp) const {
	const double	cp = Cp[i-1];
	const int	hp1 = hp % NH();
	const int	hp2 = hp / NH();
	const int	prev_hp1 = prev_hp % NH();
	const int	prev_hp2 = prev_hp / NH();
	return log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
		   log(prev_hp2 != hp2 ? cp : 1.0 - cp);
}

double ParentProgenyImputer::transition_probability(size_t i,
													int prev_h, int h) const {
	const size_t	N = num_progenies();
	const int	prev_hc = prev_h & ((1 << (N*2)) - 1);
	const int	hc = h & ((1 << (N*2)) - 1);
	const int	prev_hp = prev_h >> (N*2);
	const int	hp = h >> (N*2);
	return progeny_transition_probability(i, prev_hc, hc) +
			parent_transition_probability(i, prev_hp, hp);
}

vector<ParentProgenyImputer::DP> ParentProgenyImputer::initialize_dp() const {
	const size_t	L = NH() * NH() << (2*num_progenies());
	vector<DP>	dp(M(), DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = ref_records[0];
	const size_t	phased_index = is_mat_imputed ? 0 : 1;
	const size_t	non_phased_index = is_mat_imputed ? 1 : 0;
	const int	phased_parent_gt = record->get_geno()[phased_index] & 3;
	const int	op = record->unphased(non_phased_index);
	vector<int>	ocs;	// observed progenies' genotypes
	for(size_t ic = 0; ic != num_progenies(); ++ic) {
		ocs.push_back(record->unphased(ic+2));
	}
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability(0, h, op, ocs,
														phased_parent_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
ParentProgenyImputer::collect_possible_previous_hidden_states() const {
	const size_t	N = num_progenies();
	const size_t	L = NH() * NH() << (2*N);
	const size_t	Lc = 1 << (2*N);	// number of children states
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hc = h & (Lc - 1);
		const int	hp = h >> (N*2);
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		// Neither side will crossover
		vector<int>	prev_hps;
		for(int h1 = 0; h1 < (int)NH(); ++h1) {
			const int	prev_h1 = h1 + hp2 * NH();
			prev_hps.push_back(prev_h1);
		}
		for(int h2 = 0; h2 < (int)NH(); ++h2) {
			if(h2 == hp2)	// for duplication
				continue;
			const int	prev_h1 = hp1 + h2 * NH();
			prev_hps.push_back(prev_h1);
		}
		
		vector<int>	prev_hcs;
		for(int hc1 = 0; hc1 < 1 << (N*2); ++hc1) {
			const int	t = hc ^ hc1;	// each bit means whether crossover
			int	num_crossovers = 0;
			for(int k = 0; k < (int)N*2; ++k) {
				if(((t >> k) & 1) == 1)
					num_crossovers += 1;
			}
			if(num_crossovers <= 1)
				prev_hcs.push_back(hc1);
		}
		
		for(auto p = prev_hps.begin(); p != prev_hps.end(); ++p) {
			for(auto q = prev_hcs.begin(); q != prev_hcs.end(); ++q) {
				if(*p == hp || *q == hc) {
					const int	prev_h = (*p << (N*2)) | *q;
					prev_h_table[h].push_back(prev_h);
				}
			}
		}
	}
	
	return prev_h_table;
}

int ParentProgenyImputer::compute_non_phased_parent_gt(int h, size_t i) {
	const size_t	N = this->num_progenies();
	const int	hp = h >> (N*2);
	const int	hp1 = hp % this->NH();
	const int	hp2 = hp / this->NH();
	return this->ref_haps[hp1][i] | (this->ref_haps[hp2][i] << 1);
}

void ParentProgenyImputer::update_dp(size_t i, vector<DP>& dp) const {
	const size_t	N = num_progenies();
	const size_t	L = NH() * NH() << (2*N);
	const VCFFamilyRecord	*record = ref_records[i];
	const size_t	phased_index = is_mat_imputed ? 0 : 1;
	const size_t	non_phased_index = is_mat_imputed ? 1 : 0;
	const int	phased_parent_gt = record->get_geno()[phased_index] & 3;
	const int	op = record->unphased(non_phased_index);
	vector<int>	ocs;	// observed progenies' genotypes
	for(size_t ic = 0; ic != num_progenies(); ++ic) {
		ocs.push_back(record->unphased(ic+2));
	}
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, op, ocs,
														phased_parent_gt);
		
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			
			const double	T = transition_probability(i, prev_h, h);
			const double	prob = dp[i-1][prev_h].first + (T + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void ParentProgenyImputer::update_genotypes(const vector<int>& hs) {
	const size_t	N = num_progenies();
	for(size_t i = 0; i < this->M(); ++i) {
		VCFFamilyRecord	*record = records[i];
		int	mat_gt = 0;
		int	pat_gt = 0;
		if(this->is_mat_imputed) {
			pat_gt = this->compute_non_phased_parent_gt(hs[i], i);
			record->set_geno(1, pat_gt | 4);
			mat_gt = record->get_geno()[0] & 3;
		}
		else {
			mat_gt = this->compute_non_phased_parent_gt(hs[i], i);
			record->set_geno(0, mat_gt | 4);
			pat_gt = record->get_geno()[1] & 3;
		}
		for(size_t j = 0; j < N; ++j) {		// each progeny
			const int	hc = hs[i] >> (j * 2);
			const int	hc1 = hc & 1;
			const int	hc2 = (hc >> 1) & 1;
			const int	gtc_int = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt);
			record->set_geno(j+2, gtc_int | 4);
		}
	}
}

void ParentProgenyImputer::impute() {
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->M(); ++i) {
		update_dp(i, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}
