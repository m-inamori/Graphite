#ifndef __PARENTSIMPUTERBYPROGENY
#define __PARENTSIMPUTERBYPROGENY

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ParentsImputerByProgeny ////////////////////

class ParentsImputerByProgeny : public VCFHMM<VCFFamilyRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>		ref_haps_mat;
	const std::vector<std::vector<int>>		ref_haps_pat;
	const std::vector<std::vector<int>>		prev_h_table;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	ParentsImputerByProgeny(const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_hs_mat,
							const std::vector<std::vector<int>>& ref_hs_pat,
							const Map& map_, double w);
	~ParentsImputerByProgeny() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<VCFFamilyRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<VCFFamilyRecord *>& rs) const;
	std::size_t NH() const { return ref_haps_mat.size(); }
	std::size_t num_states() const { return (NH()*NH()) << 3; }
	std::size_t M() const { return ref_haps_mat[0].size(); }
	std::size_t num_progenies() const {
		return records[0]->num_samples() - 2;
	}
	
	int gt_by_haplotypes(int hc1, int hc2, int mat_gt, int pat_gt) const {
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
	}
	
	std::pair<int, int> compute_parent_gt(int i, int h) const;
	double emission_probability(std::size_t i, int h) const;
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	// Collect all possible previous hidden states for each hidden state.
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	int compute_non_phased_parent_gt(int h, std::size_t i);
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
