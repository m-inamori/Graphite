#include <cmath>
#include "../include/VCFOneParentImputedRough.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputedRough::VCFOneParentImputedRough(
							const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat_imputed_, const Map& map_, double w) :
				VCFBase(header, s), VCFSmallBase(),
				VCFFamilyBase(), VCFMeasurable(map_),
				records(rs), ref_haps(ref_hs), is_mat_imputed(is_mat_imputed_),
				E{{log(1.0-w*3), log(w),       log(w),       log(w)},
				  {log(w),       log(1.0-w*3), log(w),       log(w)},
				  {log(w),       log(1.0-w*3), log(w),       log(w)},
				  {log(w),       log(w),       log(1.0-w*3), log(w)}},
				Epc{{{log(1.0-w*3),  log(w),       log(w),        log(w)},
					 {log(0.5-w),    log(0.5-w),   log(w),        log(w)},
					 {log(w),        log(1.0-w*3), log(w),        log(w)}},
					{{log(0.5-w),    log(0.5-w),   log(w),        log(w)},
					 {log(0.25-w/3), log(0.5-w/3), log(0.25-w/3), log(w)},
					 {log(w),        log(0.5-w),   log(0.5-w),    log(w)}},
					{{log(w),        log(1.0-w*3), log(w),        log(w)},
					 {log(w),        log(0.5-w),   log(0.5-w),    log(w)},
					 {log(w),        log(w),       log(1.0-w*3),  log(w)}}} { }

VCFOneParentImputedRough::~VCFOneParentImputedRough() {
	Common::delete_all(records);
}

int VCFOneParentImputedRough::compute_non_phased_parent_gt(int h, int i) const {
	const size_t	NH = this->ref_haps.size();
	const int	hp1 = h % NH;
	const int	hp2 = h / NH;
	return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
}

double VCFOneParentImputedRough::emission_probability_for_parent(
												size_t i, int h, int op,
												const vector<int>& ocs,
												int phased_parent_gt) const {
	const size_t	N = num_progenies();
	const size_t	NH = this->ref_haps.size();
	const int	hp1 = h % NH;
	const int	hp2 = h / NH;
	const auto	*record = records[i];
	const int	non_phased_parent_gt = compute_non_phased_parent_gt(h, i);
	// emission probability of parent
	const double	Ep = E[non_phased_parent_gt][op];
	double	Ec = 0.0;	// emission probability of progenies
	const int	non_phased_gt = ref_haps[hp1][i] + ref_haps[hp2][i];
	const int	phased_gt = (phased_parent_gt & 1) + (phased_parent_gt >> 1);
	for(size_t j = 0; j < N; ++j) {
		const int	oc = Genotype::gt_to_int(record->get_gt(j+2));
		Ec += Epc[non_phased_gt][phased_gt][oc];
	}
	return Ep + Ec;
}

vector<VCFOneParentImputedRough::DP>
				VCFOneParentImputedRough::initialize_parent_dp() const {
	const size_t	NH = this->ref_haps.size();
	const size_t	L = NH * NH;
	const size_t	M = this->size();
	vector<DP>	dp(M, DP(L, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = records[0];
	const int	phased_parent_gt = Genotype::phased_gt_to_int(
												get_phased_parent_gt(record));
	// observed parent
	const int	op = Genotype::gt_to_int(get_non_phased_parent_gt(record));
	
	// observed progs
	vector<int>	ocs;
	for(size_t i = 0; i != num_progenies(); ++i)
		ocs.push_back(Genotype::gt_to_int(record->get_gt(i+2)));
	
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const double	E_all = emission_probability_for_parent(0, h, op, ocs,
															phased_parent_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

// Collect possible previous hidden states for the hidden state.
vector<vector<int>>
VCFOneParentImputedRough::collect_possible_previous_hidden_states() const {
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

// -> Morgan
double VCFOneParentImputedRough::dist(const VCFRecord *r1,
									  const VCFRecord *r2) const {
	const double	d = (cM(r2->pos()) - cM(r1->pos())) / 100;
	if(d != 0.0)
		return d;
	else	// probably outside map
		return (r2->pos() - r1->pos()) * 1e-6;
}

void VCFOneParentImputedRough::update_dp_for_parent(size_t i, vector<DP>& dp,
								const vector<vector<int>>& prev_h_table) const {
	const size_t	NH = this->ref_haps.size();
	const size_t	L = NH * NH;
	// Because the parent is not a direct descendant of the reference,
	// it is easier to cross over than a direct descendant.
	// its multiplication factor
	const double	K = 5.0;
	const VCFFamilyRecord	*record = records[i];
	const double	morgan = dist(records[i-1], record);
	const double	cp = Map::Kosambi(morgan * K);
	
	const int	phased_parent_gt = Genotype::phased_gt_to_int(
												get_phased_parent_gt(record));
	// observed parent
	const int	op = Genotype::gt_to_int(get_non_phased_parent_gt(record));
	// observed progenies
	vector<int>	ocs;
	for(size_t i = 0; i < num_progenies(); ++i) {
		ocs.push_back(Genotype::gt_to_int(record->get_gt(i+2)));
	}
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability_for_parent(i, h, op, ocs,
															phased_parent_gt);
		
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

vector<int> VCFOneParentImputedRough::trace_back(const vector<DP>& dps) const {
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

void VCFOneParentImputedRough::set_non_phased_parent_gt(int gt,
													VCFFamilyRecord *record) {
	const int	i = is_mat_imputed ? 1 : 0;
	const string	GT = Genotype::int_to_phased_gt(gt);
	record->set_GT(i, GT);
}

void VCFOneParentImputedRough::update_parent_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->size(); ++i) {
		const int	h = hs[i];
		VCFFamilyRecord	*record = records[i];
		const int	non_phased_parent_gt = compute_non_phased_parent_gt(h, i);
		set_non_phased_parent_gt(non_phased_parent_gt, record);
	}
}

void VCFOneParentImputedRough::impute_parent() {
	const auto	prev_h_table = collect_possible_previous_hidden_states();
	
	// DP
	vector<DP>	dp = initialize_parent_dp();
	for(int i = 1; i < (int)this->size(); ++i) {
		update_dp_for_parent(i, dp, prev_h_table);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_parent_genotypes(hs);
}

int VCFOneParentImputedRough::gt_by_haplotypes(int hc1, int hc2,
												int mat_gt, int pat_gt) {
	return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
}

double VCFOneParentImputedRough::emission_probability_for_progeny(int h, int oc,
												int mat_gt, int pat_gt) const {
	const int	hc1 = h & 1;
	const int	hc2 = h >> 1;
	const int	gtc = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt);
	return E[gtc][oc];
}

vector<VCFOneParentImputedRough::DP>
			VCFOneParentImputedRough::initialize_progeny_dp(size_t ic) const {
	const size_t	M = this->size();
	vector<DP>	dp(M, DP(4, pair<double, int>(MIN_PROB, 0)));
	const VCFFamilyRecord	*record = records[0];
	const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
	// observed progeny
	const int	oc = Genotype::gt_to_int(record->get_gt(ic+2));
	for(int h = 0; h < 4; ++h) {	// hidden state
		const double	E_all = emission_probability_for_progeny(h, oc,
																mat_gt, pat_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

void VCFOneParentImputedRough::update_dp_for_progeny(size_t i, size_t ic,
														vector<DP>& dp) const {
	const VCFFamilyRecord	*record = records[i];
	const double	morgan = dist(records[i-1], record);
	const double	cc = Map::Kosambi(morgan);
	const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
	const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
	// observed progeny
	const int	oc = Genotype::gt_to_int(record->get_gt(ic+2));
	
	for(int h = 0; h < 4; ++h) {
		const double	E_all = emission_probability_for_progeny(h, oc,
																mat_gt, pat_gt);
		
		for(int prev_h = 0; prev_h < 4; ++prev_h) {
			const int	t1 = prev_h ^ h;
			
			// log of transition probability for parent
			const double	Tc = log((t1 & 1) == 1 ? cc : 1.0 - cc) +
								 log((t1 >> 1) == 1 ? cc : 1.0 - cc);
			const double	prob = dp[i-1][prev_h].first + (Tc + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

void VCFOneParentImputedRough::update_progeny_genotypes(
										const vector<int>& hs, int ic) {
	for(size_t i = 0; i < this->size(); ++i) {
		const int	h = hs[i];
		VCFFamilyRecord	*record = records[i];
		const int	mat_gt = Genotype::phased_gt_to_int(record->get_gt(0));
		const int	pat_gt = Genotype::phased_gt_to_int(record->get_gt(1));
		const int	hc1 = h & 1;
		const int	hc2 = h >> 1;
		const int	gtc_int = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt);
		record->set_GT(ic+2, Genotype::int_to_phased_gt(gtc_int));
	}
}

void VCFOneParentImputedRough::impute_progeny(size_t ic) {
	// DP
	vector<DP>	dp = initialize_progeny_dp(ic);
	for(int i = 1; i < (int)this->size(); ++i) {
		update_dp_for_progeny(i, ic, dp);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_progeny_genotypes(hs, ic);
}

void VCFOneParentImputedRough::impute() {
	impute_parent();
	for(size_t ic = 0; ic < num_progenies(); ++ic) {
		impute_progeny(ic);
	}
}
