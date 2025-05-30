#ifndef __TYPEDETERMINER
#define __TYPEDETERMINER

#include <vector>
#include <map>
#include <tuple>


//////////////////// ParentComb ////////////////////

enum class ParentComb {
	P00x00 = 0, P00x01, P01x01, P00x11, P01x11, P11x11, PNA
};


//////////////////// TypeDeterminer ////////////////////

class TypeDeterminer {
	typedef std::tuple<std::size_t,std::size_t,std::size_t>	State;
	typedef std::tuple<double,std::size_t,std::size_t,std::size_t>	PQState;
	
	const std::size_t	N;
	const double	alpha;
	std::map<State, std::vector<std::pair<double, ParentComb>>>	memo;
	
public:
	TypeDeterminer(std::size_t n, double alpha_);
	
	std::vector<std::pair<double, ParentComb>> determine(
												const State& counter) const;
	
private:
	void make_memo00();
	void make_memo01();
	void make_memo02();
	void make_memo11();
	void make_memo12();
	void make_memo22();
	void sort();
	
	std::vector<std::pair<int, double>> binomial(std::size_t M) const;
	double genotype_probability(std::size_t n0,
								std::size_t n1, std::size_t n2) const;
	PQState initialize_state(size_t M) const;
	PQState create_pqstate(std::size_t n0,
							std::size_t n1, std::size_t n2) const;
	PQState create_pqstate(const State& s) const;
	State create_state(const PQState& s) const;
	State reverse_state(const State& s) const;
	std::vector<PQState> neighbor_states(const PQState& s0) const;
	void insert(const State& s, ParentComb c, double p);
	void insert(std::size_t n0, std::size_t n1, std::size_t n2,
													ParentComb c, double p);
	
public:
	static bool is_same_parent_gts(ParentComb c) {
		return c == ParentComb::P00x00 || c == ParentComb::P01x01 ||
											c == ParentComb::P11x11;
	}
	static bool is_homohomo(ParentComb pc) {
		return pc == ParentComb::P00x00 || pc == ParentComb::P00x11 ||
											pc == ParentComb::P11x11;
	}
	static bool is_heterohomo(ParentComb pc) {
		return pc == ParentComb::P00x01 || pc == ParentComb::P01x11;
	}
	static int get_avoiding_gt(ParentComb c) {
		return (5 - static_cast<int>(c)) >> 1;
	}
	
	static std::pair<int,int> int_gt_pair(ParentComb comb);

};
#endif
