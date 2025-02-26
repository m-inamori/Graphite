#ifndef __PARENTIMPUTER
#define __PARENTIMPUTER

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ParentImputer ////////////////////

class ParentImputer : public VCFHMM<VCFFamilyRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const bool is_mat_imputed;
	const std::vector<double>	Cp;
	const double	Epc[3][4][4];	// exhaust probability for progenies
	
public:
	ParentImputer(const std::vector<VCFFamilyRecord *>& records,
						bool is_mat_imputed,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w);
	~ParentImputer() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cp(const std::vector<VCFFamilyRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	std::size_t phased_col() const { return is_mat_imputed ? 9 : 10; }
	std::size_t num_progenies() const { return records[0]->num_samples() - 2; }
	
	int compute_phased_gt_by_refhaps(int hp, std::size_t i) const {
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	double progs_emission_probability(std::size_t i,
										int mat_gt, int pat_gt) const;
	double mat_emission_probability(std::size_t i, int h,
										int mat_gt, int pat_gt) const;
	double pat_emission_probability(std::size_t i, int h,
										int mat_gt, int pat_gt) const;
	double emission_probability(std::size_t i, int h,
										int mat_gt, int pat_gt) const;
	
	double parent_transition_probability(std::size_t i,
											int prev_hp, int hp) const;
	
	std::vector<DP> initialize_dp() const;
	
	// hidden stateに対して、可能な前のhidden stateを集めておく
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
