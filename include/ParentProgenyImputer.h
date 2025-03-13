#ifndef __PARENTPROGENYIMPUTER
#define __PARENTPROGENYIMPUTER

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ParentProgenyImputer ////////////////////

class ParentProgenyImputer : public VCFHMM<VCFFamilyRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	ref_records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const bool is_mat_imputed;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	ParentProgenyImputer(const std::vector<VCFFamilyRecord *>& rs,
						bool is_mat_imputed,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w);
	~ParentProgenyImputer() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<VCFFamilyRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<VCFFamilyRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	std::size_t num_progenies() const {
		return ref_records[0]->num_samples() - 2;
	}
	
	int gt_by_haplotypes(int hc1, int hc2, int mat_gt, int pat_gt) const {
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
	}
	
	int compute_phased_gt_by_refhaps(int hp, std::size_t i) const {
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	double emission_probability(std::size_t i, int h, int op,
								const std::vector<int>& ocs,
								int phased_parent_gt) const;
	
	double progeny_transition_probability(std::size_t i,
											int prev_hc, int hc) const;
	double parent_transition_probability(std::size_t i,
											int prev_hp, int hp) const;
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	// hidden stateに対して、可能な前のhidden stateを集めておく
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
