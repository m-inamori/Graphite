#ifndef __ORPHANIMPUTER
#define __ORPHANIMPUTER

#include "VCF.h"
#include "VCFHMM.h"
#include "Map.h"


//////////////////// OrphanImputer ////////////////////

class OrphanImputer : public VCFHMM<VCFRecord> {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFRecord *>&	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const std::vector<std::vector<int>>	prev_h_table;
	const std::vector<double>	Cp;
	
public:
	OrphanImputer(const std::vector<VCFRecord *>& records,
					const std::vector<std::vector<int>>& ref_haps,
					const Map& map_, double w);
	~OrphanImputer() { }
	
	void impute(std::size_t io);
	
private:
	std::vector<double> calc_Cp(const std::vector<VCFRecord *>& rs) const;
	std::vector<std::vector<int>>
						collect_possible_previous_hidden_states() const;
	
	std::size_t NH() const { return ref_haps.size(); }
	std::size_t M() const { return records.size(); }
	std::size_t num_progenies() const { return records[0]->num_samples() - 2; }
	
	int compute_phased_gt_by_refhaps(int hp, std::size_t i) const {
		const int	hp1 = hp % NH();
		const int	hp2 = hp / NH();
		return ref_haps[hp1][i] | (ref_haps[hp2][i] << 1);
	}
	
	double emission_probability(std::size_t i, int h, int orphan_gt) const;
	
	double transition_probability(std::size_t i, int prev_h, int h) const {
		const double	cp = Cp[i-1];
		const int	h1 = h % NH();
		const int	h2 = h / NH();
		const int	prev_h1 = prev_h % NH();
		const int	prev_h2 = prev_h / NH();
		return (prev_h1 != h1 ? log(cp) : log(1.0 - cp)) +
			   (prev_h2 != h2 ? log(cp) : log(1.0 - cp));
	}
	
	std::vector<DP> initialize_dp(std::size_t io) const;
	void update_dp(std::size_t i, std::size_t io, std::vector<DP>& dp) const;
	void update_genotypes(const std::vector<int>& hs, std::size_t io);
};
#endif
