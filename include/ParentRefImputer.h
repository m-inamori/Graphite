#ifndef __PARENTREFIMPUTER
#define __PARENTREFIMPUTER

#include <array>

#include "VCFFamily.h"
#include "VCFHMMRef.h"
#include "Map.h"


//////////////////// ParentRefImputer ////////////////////

class ParentRefImputer : public VCFHMMRef {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	records;
	const bool	is_mat_ref;
	const std::vector<std::vector<int>>	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const std::array<std::array<std::array<double, 5>, 4>, 4>	Epc;
	const std::vector<double>	Cp;
	
public:
	ParentRefImputer(const std::vector<VCFFamilyRecord *>& rs,
								const Map& map_, bool is_mat_ref_,
								const VCFGeno *ref_vcf, double w) :
						VCFHMMRef(map_, w),
						records(rs),
						is_mat_ref(is_mat_ref_),
						ref_haps(ref_vcf->create_ref_haps()),
						prev_h_table(collect_possible_previous_hidden_states()),
						Epc(calc_Epc(w)),
						Cp(calc_Cp()) { }
						
	~ParentRefImputer() { }
	
	void impute();
	
private:
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_haps[0].size(); }
	std::size_t phased_index() const { return is_mat_ref ? 1 : 0; }
	std::size_t num_progenies() const { return records[0]->num_samples()-2; }
	
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
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
	
	static std::array<std::array<std::array<double, 5>, 4>, 4>
	calc_Epc(double w);
	std::vector<double> calc_Cp() const;
};
#endif
