#ifndef __SELFPARENTIMPUTERLESSIMPUTED
#define __SELFPARENTIMPUTERLESSIMPUTED

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// SelfParentImputerLessImputed ////////////////////

class SelfParentImputerLessImputed : public VCFHMM<VCFRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFRecord *>&	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const std::vector<double>	Cp;
	const double	Epc[4][4];
	
public:
	SelfParentImputerLessImputed(const std::vector<VCFRecord *>& rs,
								 const std::vector<std::vector<int>>& ref_hs,
								 const Map& map_, double w);
	~SelfParentImputerLessImputed() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cp(const std::vector<VCFRecord *>& rs) const;
	std::vector<std::vector<int>> collect_possible_previous_hidden_states();
	
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return records.size(); }
	std::size_t num_states() const { return NH() * NH(); }
	std::size_t num_progenies() const {
		return records[0]->num_samples() - 1;
	}
	
	int parent_genotype(int hp, int i) const {
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	double progs_emission_probability(const std::vector<int>& ocs,
												int parent_gt) const;
	double emission_probability(std::size_t i, int h, int op,
										const std::vector<int>& ocs) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
