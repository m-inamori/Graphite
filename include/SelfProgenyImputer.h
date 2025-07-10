#ifndef __SELFPROGENYIMPUTER
#define __SELFPROGENYIMPUTER

#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// SelfProgenyImputer ////////////////////

class SelfProgenyImputer : public VCFHMM<VCFRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFRecord *>&	records;
	const std::size_t	ic;
	const std::vector<double>	Cc;
	
public:
	SelfProgenyImputer(const std::vector<VCFRecord *>& rs,
						std::size_t iprog, const Map& map_, double w);
	~SelfProgenyImputer() { }
	
	void impute();
	
private:
	std::vector<double> calc_Cc(const std::vector<VCFRecord *>& rs) const;
	std::size_t M() const { return records.size(); }
	std::size_t num_progenies() const {
		return records[0]->num_samples() - 1;
	}
	
	int progeny_genotype(int h, int parent_gt) const {
		const int	hc1 = h & 1;
		const int	hc2 = h >> 1;
		return ((parent_gt >> hc1) & 1) | (((parent_gt >> hc2) & 1) << 1);
	}
	
	double emission_probability(int h, int parent_gt, int oc) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const;
	
	std::vector<DP> initialize_dp() const;
	
	void update_dp(std::size_t i, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs);
};
#endif
