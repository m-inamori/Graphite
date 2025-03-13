#include <cmath>
#include "../include/VCFOneParentImputed.h"
#include "../include/common.h"

using namespace std;

VCFOneParentImputed::VCFOneParentImputed(const std::vector<STRVEC>& header,
							const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs,
							bool is_mat_imputed_, const Map& map_, double w) :
				VCFOneParentImputedBase(header, s, rs), VCFMeasurable(map_),
				ref_haps(ref_hs), is_mat_imputed(is_mat_imputed_),
				E{{log(1.0-w*2), log(w/2),     log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(w/2),     log(1.0-w*2), log(w)}} { }

VCFOneParentImputed::~VCFOneParentImputed() {
	Common::delete_all(records);
}

int VCFOneParentImputed::gt_by_haplotypes(int hc1, int hc2,
											int phased_parent_gt,
											int non_phased_parent_gt) const {
	const int	mat_gt = is_mat_imputed ? phased_parent_gt
										: non_phased_parent_gt;
	const int	pat_gt = is_mat_imputed ? non_phased_parent_gt
										: phased_parent_gt;
	return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
}

int VCFOneParentImputed::compute_non_phased_parent_gt(int h, int i) const {
	const size_t	N = this->num_progenies();
	const size_t	NH = this->ref_haps.size();
	const int	hp = h >> (N*2);
	const int	hp1 = hp % NH;
	const int	hp2 = hp / NH;
	return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
}

double VCFOneParentImputed::emission_probability(size_t i, int h, int op,
												 const vector<int>& ocs,
												 int phased_parent_gt) const {
	const size_t	N = num_progenies();
	const int	non_phased_parent_gt = compute_non_phased_parent_gt(h, i);
	// emission probability of parent
	const double	Ep = E[non_phased_parent_gt][op];
	double	Ec = 0.0;	// emission probability of progenies
	for(size_t j = 0; j < N; ++j) {
		const int	hc1 = (h >> (j * 2)) & 1;
		const int	hc2 = (h >> (j * 2 + 1)) & 1;
		const int	gtc = gt_by_haplotypes(hc1, hc2, phased_parent_gt,
														non_phased_parent_gt);
		Ec += E[gtc][ocs[j]];
	}
	return Ep + Ec;
}

vector<VCFOneParentImputed::DP> VCFOneParentImputed::initialize_dp() const {
	const size_t	L = compute_num_hidden_states();
	const size_t	M = this->size();
	vector<DP>	dp(M, DP(L, pair<double, int>(-1e300, 0)));
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
		const double	E_all = emission_probability(0, h, op, ocs,
														phased_parent_gt);
		dp[0][h] = make_pair(E_all, h);
	}
	return dp;
}

bool VCFOneParentImputed::is_few_crossover(int t) {
	if(t == 0)	// no crossover
		return true;
	
	int	i = 0;
	while(((t >> i) & 1) == 0) {
		++i;
	}
	return t == 1 << i;
}

// hidden stateに対して、可能な前のhidden stateを集めておく
vector<vector<int>>
		VCFOneParentImputed::collect_possible_previous_hidden_states() const {
	const size_t	L = compute_num_hidden_states();
	const size_t	N = this->num_progenies();
	const size_t	Lc = 1 << (2*N);
	const size_t	NH = this->ref_haps.size();
	vector<vector<int>>	prev_h_table(L, vector<int>());
	for(int h = 0; h < (int)L; ++h) {	// hidden state
		const int	hc = h & (Lc - 1);
		const int	hp = h >> (N*2);
		const int	hp1 = hp % NH;
		const int	hp2 = hp / NH;
		// 両側乗り越えることはないとする
		// non-phasedの親のあり得る前の状態
		vector<int>	prev_hps;
		for(int h1 = 0; h1 < (int)NH; ++h1) {
			const int	prev_h1 = h1 + hp2 * NH;
			prev_hps.push_back(prev_h1);
		}
		for(int h2 = 0; h2 < (int)NH; ++h2) {
			if(h2 == hp2)	// for duplication
				continue;
			const int	prev_h1 = hp1 + h2 * NH;
			prev_hps.push_back(prev_h1);
		}
		
		// 後代のあり得る前の状態
		vector<int>	prev_hcs;
		for(int hc1 = 0; hc1 < 1 << (N*2); ++hc1) {
			const int	t = hc ^ hc1;	// 乗り換えたかをビットで表す
			// 乗り換えは合わせて1回まで許容できる
			if(is_few_crossover(t)) {
				prev_hcs.push_back(hc1);
			}
		}
		
		for(auto p = prev_hps.begin(); p != prev_hps.end(); ++p) {
			const int	prev_hp = *p;
			for(auto q = prev_hcs.begin(); q != prev_hcs.end(); ++q) {
				const int	prev_hc = *q;
				// 両方乗り換えているのはダメ
				if(prev_hp == hp || prev_hc == hc) {
					const int	prev_h = (prev_hp << (N*2)) | prev_hc;
					prev_h_table[h].push_back(prev_h);
				}
			}
		}
	}
	
	return prev_h_table;
}

// -> Morgan
double VCFOneParentImputed::dist(const VCFRecord *r1,
								 const VCFRecord *r2) const {
	const double	d = (cM(r2->pos()) - cM(r1->pos())) / 100;
	if(d != 0.0)
		return d;
	else	// probably outside map
		return (r2->pos() - r1->pos()) * 1e-6;
}

void VCFOneParentImputed::update_dp(size_t i, vector<DP>& dp,
								const vector<vector<int>>& prev_h_table) const {
	const size_t	N = this->num_progenies();
	const size_t	NH = this->ref_haps.size();
	const size_t	L = compute_num_hidden_states();
	// Because the parent is not a direct descendant of the reference,
	// it is easier to cross over than a direct descendant.
	// its multiplication factor
	const double	K = 5.0;
	const VCFFamilyRecord	*record = records[i];
	const double	morgan = dist(records[i-1], record);
	const double	cc = Map::Kosambi(morgan);
	const double	cp = Map::Kosambi(morgan * K);
	
	const int	phased_parent_gt = Genotype::phased_gt_to_int(
												get_phased_parent_gt(record));
	// observed parent
	const int	op = Genotype::gt_to_int(get_non_phased_parent_gt(record));
	// observed progenies
	vector<int>	ocs;
	for(size_t ic = 0; ic < num_progenies(); ++ic) {
		ocs.push_back(Genotype::gt_to_int(record->get_gt(ic+2)));
	}
	
	for(int h = 0; h < (int)L; ++h) {
		const double	E_all = emission_probability(i, h, op, ocs,
														phased_parent_gt);
		
		const int	hp = h >> (N*2);
		const int	hp1 = hp % NH;
		const int	hp2 = hp / NH;
		const auto&	prev_hs = prev_h_table[h];
		for(auto p = prev_hs.begin(); p != prev_hs.end(); ++p) {
			const int	prev_h = *p;
			const int	prev_hp = prev_h >> (N*2);
			const int	t = prev_h ^ h;		// whether crossover or not in bits
			
			// log of transition probability for progenies
			double	Tc = 0.0;
			for(size_t j = 0; j < N*2; ++j) {
				Tc += log(((t >> j) & 1) == 1 ? cc : 1.0 - cc);
			}
			
			// log of transition probability for parent
			const int	prev_hp1 = prev_hp % NH;
			const int	prev_hp2 = prev_hp / NH;
			const double	Tp = log(prev_hp1 != hp1 ? cp : 1.0 - cp) +
								 log(prev_hp2 != hp2 ? cp : 1.0 - cp);
			const double	prob = dp[i-1][prev_h].first + (Tc + Tp + E_all);
			dp[i][h] = max(dp[i][h], make_pair(prob, prev_h));
		}
	}
}

vector<int> VCFOneParentImputed::trace_back(const vector<DP>& dps) const {
	const size_t	M = this->size();
	vector<int>	hs(M, 0);
	
	pair<double, int>	max_pair(-1e300, 0);
	int	max_h = 0;
	const size_t	L = compute_num_hidden_states();
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

void VCFOneParentImputed::set_non_phased_parent_gt(int gt,
													VCFFamilyRecord *record) {
	const int	i = is_mat_imputed ? 1 : 0;
	const string	GT = Genotype::int_to_phased_gt(gt);
	record->set_GT(i, GT);
}

void VCFOneParentImputed::update_genotypes(const vector<int>& hs) {
	for(size_t i = 0; i < this->size(); ++i) {
		const int	h = hs[i];
		VCFFamilyRecord	*record = records[i];
		const int	phased_parent_gt = Genotype::phased_gt_to_int(
												get_phased_parent_gt(record));
		const int	non_phased_parent_gt = compute_non_phased_parent_gt(h, i);
		set_non_phased_parent_gt(non_phased_parent_gt, record);
		
		for(size_t j = 0; j < num_progenies(); ++j) {
			const int	hc1 = (hs[i] >> (j * 2)) & 1;
			const int	hc2 = (hs[i] >> (j * 2 + 1)) & 1;
			const int	gtc = gt_by_haplotypes(hc1, hc2, phased_parent_gt,
														non_phased_parent_gt);
			record->set_GT(j+2, Genotype::int_to_phased_gt(gtc));
		}
	}
}

void VCFOneParentImputed::impute() {
	const auto	prev_h_table = collect_possible_previous_hidden_states();
	
	// DP
	vector<DP>	dp = initialize_dp();
	for(int i = 1; i < (int)this->size(); ++i) {
		update_dp(i, dp, prev_h_table);
	}
	
	const vector<int>	hs = trace_back(dp);
	update_genotypes(hs);
}

size_t VCFOneParentImputed::amount() const {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH*NH * (2*NH - 1);
	return R * M;
}
