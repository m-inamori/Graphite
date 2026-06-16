#ifndef __IMPUTERBYPARENTPROGENY
#define __IMPUTERBYPARENTPROGENY

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ImputerByParentProgeny ////////////////////

class ImputerByParentProgeny : public VCFHMM<VCFFamilyRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const bool should_impute_mat;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	ImputerByParentProgeny(const std::vector<VCFFamilyRecord *>& rs,
							const std::vector<std::vector<int>>& ref_haps,
							bool should_impute_mat_, const Map& map_, double w);
	~ImputerByParentProgeny() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<VCFFamilyRecord *>& rs) const;
	std::vector<double> calc_Cp(const std::vector<VCFFamilyRecord *>& rs) const;
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	std::size_t num_should_impute_progenies() const {
		return records[0]->get_genos().size() - 3;
	}
	std::size_t num_states() const {
		return NH() << (num_should_impute_progenies()*2+3);
	}
	std::size_t should_imputed_index() const {
		return static_cast<std::size_t>(!should_impute_mat);
	}
	
	void decode(int h, int& hw, int& hp1,
						int& hp21, int& hp22, std::vector<int>& hc) const;
	int parent_gt_by_haplotypes(std::size_t i, int hw,
									int hp1, int hp2, int prog_gt) const;
	int progeny_gt_by_haplotypes(std::size_t j, const std::vector<int>& hc,
												int mat_gt, int pat_gt) const;
	std::vector<int> compute_gts(std::size_t i, int h) const;
	double emission_probability(std::size_t i, int h) const;
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
