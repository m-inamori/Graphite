#ifndef __HAPLOTYPE
#define __HAPLOTYPE

#include <vector>


//////////////////// Haplotype ////////////////////

class Haplotype {
public:
	using Pair = std::pair<Haplotype, Haplotype>;
	
	const std::vector<int>	hap;
	const std::pair<std::size_t, int>	position;
	
public:
	Haplotype(const std::vector<int>& h, std::size_t sample_id, int i) :
											hap(h), position(sample_id, i) { }
	
public:
	static Haplotype default_value() {
		const std::vector<int>	hap1;
		return Haplotype(hap1, 0, 2);
	}
	
	static int score(const Pair& hap, const std::vector<int>& int_gts);
	static int score(const Haplotype& hap_mat,
					const Haplotype& hap_pat, const std::vector<int>& int_gts);
	
	static std::vector<Pair> collect_optimal_haplotype_pairs(
										const std::vector<Haplotype>& haps_mat,
										const std::vector<Haplotype>& haps_pat,
										const std::vector<int>& int_gts);
	
	static int match_score(const Pair& prev_hap, const Pair& hap);
	
	static std::vector<Pair> collect_max_score(const std::vector<Pair>& combs,
														const Pair& prev_hap);
	
	static Pair impute(const std::vector<int>& int_gts,
									const std::vector<Haplotype>& haps_mat,
									const std::vector<Haplotype>& haps_pat,
									const Pair& prev_hap, int seed);
};

using HaplotypePair = std::pair<Haplotype, Haplotype>;
#endif
