#include <vector>
#include <map>
#include <tuple>

class TypeDeterminer {
	typedef std::tuple<int,int,int>	State;
	typedef std::tuple<double,int,int,int>	PQState;
	
	const int	N;
	const double	alpha;
	// どのGenotypeの組み合わせかを4進で表す
	// e.g. 0/1 x 1/1 -> 1 + 2*4 = 9
	std::map<State,int>	memo;
	
public:
	TypeDeterminer(int n, double alpha_);
	
	int determine(const State& counter) const;
	
private:
	void make_memo00();
	void make_memo01();
	void make_memo02();
	void make_memo11();
	void make_memo12();
	void make_memo22();
	
	int binomial(int M) const;
	double genotype_probability(int n0, int n1, int n2) const;
	PQState initialize_state(int M) const;
	PQState create_pqstate(int n0, int n1, int n2) const;
	PQState create_pqstate(const State& s) const;
	State create_state(const PQState& s) const;
	State reverse_state(const State& s) const;
	std::vector<PQState> neighbor_states(const PQState& s0) const;
	void insert(const State& s, int value);
	void insert(int n0, int n1, int n2, int value);
	
public:
	static std::vector<int> int_gt_pairs(int p);
	static std::pair<int,int> int_gt_pair(int p);
};
