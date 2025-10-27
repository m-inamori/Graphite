#ifndef __SELFPARENTIMPUTER
#define __SELFPARENTIMPUTER

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// SelfParentImputer ////////////////////

class SelfParentImputer : public VCFHMM<GenoRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<GenoRecord *>&	ref_records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::size_t	ic;
	const std::vector<std::vector<int>>	prev_h_table;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	const double	Epc[4][4];	// exhaust probability for progenies
	
public:
	SelfParentImputer(const std::vector<GenoRecord *>& rs,
					  const std::vector<std::vector<int>>& ref_haps,
					  std::size_t iprog, const Map& map_, double w);
	~SelfParentImputer() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<GenoRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<GenoRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	std::size_t num_progenies() const {
		return ref_records[0]->num_samples() - 1;
	}
	
	std::tuple<int, int, int, int> decode_state(int h) const;
	std::size_t num_states() const { return (NH() + 1) * 16; }
	
	int allele(int h, std::size_t i, int prog_gt) const {
		return h < 2 ? ((prog_gt >> h) & 1) : ref_haps[h-2][i];
	}
	
	int parent_genotype(int h, std::size_t i, int prog_gt) const {
		const auto	t = decode_state(h);
		const int	hp1 = std::get<0>(t);
		const int	hp2 = std::get<1>(t);
		const int	a1 = allele(hp1, i, prog_gt);
		const int	a2 = allele(hp2, i, prog_gt);
		return a1 | (a2 << 1);
	}
	
	int progeny_genotype(int h, std::size_t i, int prog_gt) const {
		const auto	t = decode_state(h);
		const int	hc1 = std::get<2>(t);
		const int	hc2 = std::get<3>(t);
		const int	gt_parent = parent_genotype(h, i, prog_gt);
		return ((gt_parent >> hc1) & 1) | ((gt_parent >> hc2) & 1);
	}
	
	// Probability of phased genotypes being exhausted for phased genotypes
	double phased_emission_probability(int h, std::size_t i, int prog_gt) const;
	double non_phased_emission_probability(const std::vector<int>& ocs,
														int parent_gt) const;
	
	int compute_phased_gt_by_refhaps(int hp, std::size_t i) const {
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	double emission_probability(int h, std::size_t i, int op, int prog_gt,
											const std::vector<int>& ocs) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	// hidden stateに対して、可能な前のhidden stateを集めておく
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
