#ifndef __PROGENYSELFREFIMPUTER
#define __PROGENYSELFREFIMPUTER

#include <cmath>
#include "GenoRecord.h"
#include "VCFHMMSelfRef.h"
#include "Map.h"


//////////////////// ProgenySelfRefImputer ////////////////////

class ProgenySelfRefImputer : public VCFHMMSelfRef {
public:
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<GenoRecord *>&	records;
	const std::vector<double>	Cc;
	
public:
	ProgenySelfRefImputer(const std::vector<GenoRecord *>& rs,
											const Map& map_, double w);
	~ProgenySelfRefImputer() { }
	
	void impute(std::size_t j);
	
private:
	std::vector<double> calc_Cc(const std::vector<GenoRecord *>& rs) const;
	
	std::size_t M() const { return records.size(); }
	
	int gt_by_haplotypes(int hc, int parent_gt) const {
		const int	hc1 = hc & 1;
		const int	hc2 = hc >> 1;
		return ((parent_gt >> hc1) & 1) | (((parent_gt >> hc2) & 1) << 1);
	}
	
	double emission_probability(std::size_t i, std::size_t j,
											int h, int parent_gt) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const {
		const double	cc = Cc[i-1];
		const int	h1 = h & 1;
		const int	h2 = h >> 1;
		const int	prev_h1 = prev_h & 1;
		const int	prev_h2 = prev_h >> 1;
		return (prev_h1 != h1 ? log(cc) : log(1.0 - cc)) +
			   (prev_h2 != h2 ? log(cc) : log(1.0 - cc));
	}
	
	std::vector<VCFHMMSelfRef::DP> initialize_dp(std::size_t j) const;
	void update_dp(std::size_t i, std::size_t j, std::vector<DP>& dp) const;
	void update_genotypes(std::size_t j, const std::vector<int>& hs);
};
#endif
