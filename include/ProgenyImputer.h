#ifndef __PROGENYIMPUTER
#define __PROGENYIMPUTER

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ProgenyImputer ////////////////////

class ProgenyImputer : public VCFHMM {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	records;
	const std::vector<double>	Cc;
	
public:
	ProgenyImputer(const std::vector<VCFFamilyRecord *>& records,
											const Map& map_, double w);
	~ProgenyImputer() { }
	
	void impute(std::size_t j);
	
private:
	std::vector<double> calc_Cc(
						const std::vector<VCFFamilyRecord *>& rs) const;
	
	std::size_t M() const { return records.size(); }
	std::size_t num_progenies() const { return records[0]->num_samples() - 2; }
	
	int gt_by_haplotypes(int hc, int mat_gt, int pat_gt) const {
		const int	hc1 = hc & 1;
		const int	hc2 = hc >> 1;
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
	}
	
	double emission_probability(std::size_t i, std::size_t j, int h,
												int mat_gt, int pat_gt) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const {
		const double	cc = Cc[i-1];
		const int	h1 = h & 1;
		const int	h2 = h >> 1;
		const int	prev_h1 = prev_h & 1;
		const int	prev_h2 = prev_h >> 1;
		return (prev_h1 != h1 ? log(cc) : log(1.0 - cc)) +
			   (prev_h2 != h2 ? log(cc) : log(1.0 - cc));
	}
	
	std::vector<DP> initialize_dp(std::size_t j) const;
	void update_dp(std::size_t i, std::size_t j, std::vector<DP>& dp) const;
	void update_genotypes(std::size_t j, const std::vector<int>& hs);
};
#endif
