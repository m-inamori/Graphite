#ifndef __VCFHETEROIMPHOMO
#define __VCFHETEROIMPHOMO

#include "VCFHeteroHomoOnePhased.h"

class Map;


//////////////////// VCFHeteroImpHomo ////////////////////

class VCFHeteroImpHomo : public VCFHeteroHomoOnePhased {
	static const int	INF = 1000000000;
	
	////////// State //////////
	
	struct State {
		int	s;	// haplotypes & num_crossovers
		int	n;	// num of pogenies
		
		State(int s_, int n_) : s(s_), n(n_) { }
		
		int haplotype(int i) const { return (s >> i) & 1; }
		int num_crossovers() const { return s >> n; }
		bool is_full_crossovers() const { return num_crossovers() >= n*2; }
		
		void set_haplotype(int h, int i) { s |= h << i; }
		void increment_num_crossovers() { s += 1 << n; }
		void set_num_crossovers(int nc) { s |= nc << n; }
	};
	
	////////// Value //////////
	
	struct Value {
		int		num_mis;
		State	state;
		int		order;
		
		Value(int nm, State s, int o) : num_mis(nm), state(s), order(o) { }
		
		bool is_valid() const { return num_mis != INF; }
		Value add(const Value& other, State state) const {
			return Value(num_mis + other.num_mis, state, other.order);
		}
		
		bool operator <(const Value& other) const {
			if(num_mis != other.num_mis)
				return num_mis < other.num_mis;
			else if(state.num_crossovers() != other.state.num_crossovers())
				return state.num_crossovers() < other.state.num_crossovers();
			else
				return order < other.order;
		}
		
		static Value create_default(int n) {
			return Value(INF, State(0, n), 0);
		}
	};
	
	using DP = std::vector<Value>;
	
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
						const State& state, const VCFRecord *record) const;
	std::vector<std::pair<State, Value>> next_states(State state,
												const VCFRecord *record) const;
	DP update_dp(const DP& dp, const VCFRecord *record) const;
	void trace_back(State state, const std::vector<DP>& dps);
	
	void impute();
};
#endif
