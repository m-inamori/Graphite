#include "../include/Haplotype.h"

using namespace std;

int Haplotype::score(const HaplotypePair& hap, const vector<int>& int_gts) {
	int	total_score = 0;
	for(size_t k = 0; k < int_gts.size(); ++k) {
		if(hap.first.hap[k] + hap.second.hap[k] == int_gts[k])
			total_score += 1;
	}
	return total_score;
}

int Haplotype::score(const Haplotype& hap_mat,
					const Haplotype& hap_pat, const vector<int>& int_gts) {
	int	total_score = 0;
	for(size_t k = 0; k < int_gts.size(); ++k) {
		if(hap_mat.hap[k] + hap_pat.hap[k] == int_gts[k])
			total_score += 1;
	}
	return total_score;
}

int Haplotype::match_score(const HaplotypePair& prev_hap,
										const HaplotypePair& hap) {
	const Haplotype&	prev_hap1 = prev_hap.first;
	const Haplotype&	prev_hap2 = prev_hap.second;
	const Haplotype&	hap1 = hap.first;
	const Haplotype&	hap2 = hap.second;
	return (prev_hap1.position == hap1.position ? 1 : 0) +
		   (prev_hap2.position == hap2.position ? 1 : 0);
}

vector<HaplotypePair> Haplotype::collect_optimal_haplotype_pairs(
											const vector<Haplotype>& haps_mat,
											const vector<Haplotype>& haps_pat,
											const vector<int>& int_gts) {
	int	max_score = 0;
	vector<HaplotypePair>	max_combs;
	for(auto p = haps_mat.begin(); p != haps_mat.end(); ++p) {
		const Haplotype&	hap_mat = *p;
		for(auto q = haps_pat.begin(); q != haps_pat.end(); ++q) {
			const Haplotype&	hap_pat = *q;
			const int	s = score(hap_mat, hap_pat, int_gts);
			if(s > max_score) {
				max_score = s;
				max_combs.clear();
				max_combs.push_back(HaplotypePair(hap_mat, hap_pat));
			}
			else if(s == max_score) {
				max_combs.push_back(HaplotypePair(hap_mat, hap_pat));
			}
		}
	}
	return max_combs;
}

vector<HaplotypePair> Haplotype::collect_max_score(
										const vector<HaplotypePair>& combs,
										const HaplotypePair& prev_hap) {
	int	max_score = 0;
	vector<HaplotypePair>	max_combs;
	for(auto p = combs.begin(); p != combs.end(); ++p) {
		const int	score = match_score(prev_hap, *p);
		if(score == max_score) {
			max_combs.push_back(*p);
		}
		else if(score > max_score) {
			max_score = score;
			max_combs.clear();
			max_combs.push_back(*p);
		}
	}
	return max_combs;
}

HaplotypePair Haplotype::impute(const vector<int>& int_gts,
								const vector<Haplotype>& haps_mat,
								const vector<Haplotype>& haps_pat,
								const HaplotypePair& prev_hap, int seed) {
	// blute-force
	const auto	combs = collect_optimal_haplotype_pairs(haps_mat,
														haps_pat, int_gts);
	
	// collect combinations that are most matched with the previous
	const auto	filtered_combs = collect_max_score(combs, prev_hap);
	
	// select a combination with random
	const int	j = seed % filtered_combs.size();
	return filtered_combs[j];
}
