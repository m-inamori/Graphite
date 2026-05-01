#ifndef __PROGENYREFIMPUTER
#define __PROGENYREFIMPUTER

#include "VCFFamily.h"
#include "VCFHMMRef.h"

class Map;


//////////////////// ProgenyRefImputer ////////////////////

class ProgenyRefImputer : public VCFHMMRef {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	records;
	const std::vector<double>	Cc;
	
public:
	ProgenyRefImputer(const std::vector<VCFFamilyRecord *>& rs,
						const Map& map_, double w) :
										VCFHMMRef(map_, w),
										records(rs),
										Cc(calc_Cc(records)) { }
	~ProgenyRefImputer() { }
	
	void impute(std::size_t i);
	
private:
	std::size_t num_progenies() const { return records[0]->num_progenies(); }
	int gt_by_haplotypes(int hc, int mat_gt, int pat_gt) const;
	std::vector<double> calc_Cc(const std::vector<VCFFamilyRecord *>& rs) const;
	std::size_t M() const { return records.size(); }
	
	double emission_probability(std::size_t i, std::size_t j, int h,
												int mat_gt, int pat_gt) const;
	double progeny_transition_probability(std::size_t i,
											int prev_h, int h) const;
	
	std::vector<DP> initialize_dp(std::size_t j) const;
	
	void update_dp(std::size_t i, std::size_t j, std::vector<DP>& dp) const;
	void update_genotypes(std::size_t j, const std::vector<int>& hs);
	
};

#endif
