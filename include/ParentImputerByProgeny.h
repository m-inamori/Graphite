#ifndef __PARENTIMPUTERBYPROGENY
#define __PARENTIMPUTERBYPROGENY

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ParentImputerByProgeny ////////////////////

class ParentImputerByProgeny : public VCFHMM<VCFRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFRecord *>&	ref_records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const bool is_mat_known;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	ParentImputerByProgeny(const std::vector<VCFRecord *>& rs,
							const std::vector<std::vector<int>>& ref_haps,
							bool is_mat_known_, const Map& map_, double w);
	~ParentImputerByProgeny() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<VCFRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<VCFRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	
	int compute_parent_gt(std::size_t i, int h) const;
	
	double emission_probability(std::size_t i, int h) const;
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	// hidden stateに対して、可能な前のhidden stateを集めておく
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
