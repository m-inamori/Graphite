#ifndef __SELFIMPUTER
#define __SELFIMPUTER

#include <array>
#include "VCFHMM.h"
#include "GenoRecord.h"
#include "Map.h"


//////////////////// SelfImputer ////////////////////

class SelfImputer : public VCFHMM<GenoRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<GenoRecord *>&	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	SelfImputer(const std::vector<GenoRecord *>& rs,
						const std::vector<std::vector<int>>& ref_hs,
						const Map& map_, double w);
	~SelfImputer() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<GenoRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<GenoRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return records.size(); }
	std::size_t num_states() const { return NH()*NH() << (num_progenies()*2); }
	std::size_t num_progenies() const {
		return records[0]->num_samples() - 1;
	}
	std::array<int, 3> decode_state(int h) const;
	
	bool is_only_one_or_zero_crossover(int hp1, int hp2, int hc,
														int prev_h) const;
	std::vector<std::vector<int>>
	collect_possible_previous_hidden_states() const;
	
	int compute_parent_phased_gt(int h, int i) const {
		const auto	t = decode_state(h);
		const int	hp1 = std::get<0>(t);
		const int	hp2 = std::get<1>(t);
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	int parent_genotype(int hp1, int hp2, int i) const {
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	std::vector<int> compute_progeny_phased_gts(int hc, int parent_gt) const;
	double progs_emission_probability(int hc, const std::vector<int>& ocs,
														int parent_gt) const;
	double emission_probability(int i, int h, int op,
								std::vector<int>& ocs) const;
	
	double parent_transition_probability(int i, int prev_hp1, int prev_hp2,
														int hp1, int hp2) const;
	double progeny_transition_probability(int i, int prev_hc, int hc) const;
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
