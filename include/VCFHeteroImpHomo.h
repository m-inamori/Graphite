#ifndef __VCFHETEROIMPHOMO
#define __VCFHETEROIMPHOMO

#include "VCFHeteroHomoOnePhased.h"

class Map;


//////////////////// VCFHeteroImpHomo ////////////////////

class VCFHeteroImpHomo : public VCFHeteroHomoOnePhased {
	using DPValue = std::pair<int, int>;
	using DP = std::vector<DPValue>;
	const int	INF = 1000000000;
	
public:
	VCFHeteroImpHomo(const std::vector<STRVEC>& h, const STRVEC& s,
						std::vector<VCFFillableRecord *> rs,
						bool is_mat_hetero, const Map& m) :
				VCFHeteroHomoOnePhased(h, s, rs, is_mat_hetero, m) { }
	~VCFHeteroImpHomo() { }
	
	///// non-virtual methods /////
	std::size_t imputed_index() const { return is_mat_hetero ? 1 : 0; }
	std::size_t non_imputed_index() const { return is_mat_hetero ? 0 : 1; }
	
	///// DP /////
	DP init_dp() const;
	static int get_order(int state) { return state & 1; }
	static int get_hap(int state, int i) { return (state >> (i+1)) & 1; }
	int get_crossover(int state) const { return state >> (num_progenies()+1); }
	int gt_by_haplotype(int h, int gt_imputed, int order) const {
		return gt_imputed / 2 + (order ^ h);
	}
	
	std::vector<bool> is_right_gt(int order, int gt_imputed,
							int state, const VCFRecord *record) const;
	std::vector<std::pair<int, int>> next_states_small(int state,
												const VCFRecord *record) const;
	std::vector<std::pair<int, int>> next_states_large(int state,
												const VCFRecord *record) const;
	std::vector<std::pair<int, int>> next_states(int state,
												const VCFRecord *record) const;
	DP update_dp(const DP& dp, const VCFRecord *record) const;
	void trace_back(int state, const std::vector<DP>& dps);
	
	void impute();
};
#endif
