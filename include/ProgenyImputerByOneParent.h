#ifndef __PROGENYIMPUTERBYONEPARENT
#define __PROGENYIMPUTERBYONEPARENT

#include <cmath>
#include "VCFFamily.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// ProgenyImputerByOneParent ////////////////////

class ProgenyImputerByOneParent : public VCFHMM<VCFFamilyRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>&	ref_records;
	const std::vector<std::vector<int>>&	ref_haps;
	const bool	is_mat_imputed;
	const std::vector<double>	Cc;
	const std::vector<double>	Cp;
	
public:
	ProgenyImputerByOneParent(const std::vector<VCFFamilyRecord *>& rs,
								const std::vector<std::vector<int>>& ref_haps,
								bool is_mat_imputed, const Map& map_, double w);
	~ProgenyImputerByOneParent() { }
	
	void impute(std::size_t j);
	
private:
	std::vector<double> calc_Cc(
						const std::vector<VCFFamilyRecord *>& rs) const;
	std::vector<double> calc_Cp(
						const std::vector<VCFFamilyRecord *>& rs) const;
	
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return ref_records.size(); }
	std::size_t num_progenies() const {
		return ref_records[0]->num_samples() - 2;
	}
	
	std::pair<int, int> parent_alleles(int h, size_t i) const {
		const int	hc1 = h & 1;
		const int	hc2 = h >> 1;
		const int	a2 = ref_haps[hc2][i];
		if(is_mat_imputed) {
			const int	a1 = records[i]->get_mat_allele(hc1);
			return std::make_pair(a1, a2);
		}
		else {
			const int	a1 = records[i]->get_pat_allele(hc1);
			return std::make_pair(a2, a1);
		}
	}
	
	double emission_probability(std::size_t i, std::size_t j, int h,
														int parent_gt) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const {
		const double	cc = Cc[i-1];
		const double	cp = Cp[i-1];
		const int	h1 = h & 1;
		const int	h2 = h >> 1;
		const int	prev_h1 = prev_h & 1;
		const int	prev_h2 = prev_h >> 1;
		return (prev_h1 != h1 ? log(cc) : log(1.0 - cc)) +
			   (prev_h2 != h2 ? log(cp) : log(1.0 - cp));
	}
	
	std::vector<DP> initialize_dp(std::size_t j) const;
	void update_dp(std::size_t i, std::size_t j, std::vector<DP>& dp) const;
	void update_genotypes(std::size_t j, const std::vector<int>& hs);
};
#endif
