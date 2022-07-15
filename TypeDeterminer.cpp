#include <set>
#include <queue>
#include <cmath>
#include "TypeDeterminer.h"

using namespace std;


TypeDeterminer::TypeDeterminer(int n, double alpha_) : N(n), alpha(alpha_) {
	make_memo00();
	make_memo01();
	make_memo02();
	make_memo11();
	make_memo12();
	make_memo22();
}

int TypeDeterminer::binomial(int M) const {
	double	total_p = 0.0;
	double	p = pow(0.5, M);
	int	n = 1;
	for( ; total_p < alpha / 2; ++n) {
		p = p * (M - n + 1) / n;
		total_p += p;
	}
	return n;
}

void TypeDeterminer::insert(const State& s, int value) {
	if(memo.find(s) != memo.end())
		memo[s] = OVERLAP;
	else
		memo[s] = value;
}

void TypeDeterminer::insert(int n0, int n1, int n2, int value) {
	insert(State(n0, n1, n2), value);
}

void TypeDeterminer::make_memo00() {
	// あり得ないGenotypeも1割だけあってもよいことにする
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		for(int n1 = 0; n1 <= N / 10; ++n1) {
			for(int n2 = 0; n2 <= N / 10; ++n2)
				insert(N - n1 - n2 - num_NA, n1, n2, 0);
		}
	}
}

void TypeDeterminer::make_memo01() {
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		for(int n2 = 0; n2 <= N / 10; ++n2) {
			const int	M = N - num_NA - n2;
			const int	n = binomial(M);
			for(int n0 = n; n0 <= M - n; ++n0)
				insert(n0, M-n0, n2, 1);
		}
	}
}

void TypeDeterminer::make_memo02() {
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		for(int n0 = 0; n0 <= N / 10; ++n0) {
			for(int n2 = 0; n2 <= N / 10; ++n2)
				insert(n0, N - n0 - n2 - num_NA, n2, 2);
		}
	}
}

double TypeDeterminer::genotype_probability(int n0, int n1, int n2) const {
	// nが大きければ対数計算しなければならないかも
	const int	n = n0 + n1 + n2;
	const double	p0 = 0.25;
	const double	p1 = 0.5;
	const double	p2 = 0.25;
	double	p = 1.0;
	for(int i = 1; i <= n0; ++i)
		p *= p0 * (n - i + 1) / i;
	for(int i = 1; i <= n1; ++i)
		p *= p1 * (n - n0 - i + 1) / i;
	for(int i = 1; i <= n2; ++i)
		p *= p2 * (n2 - i + 1) / i;
	return p;
}

TypeDeterminer::PQState TypeDeterminer::initialize_state(int M) const {
	const int	n1 = M / 2;
	const int	n0 = (M - n1 + 1) / 2;	// n0 >= n2とする
	const int	n2 = M - n0 - n1;
	return create_pqstate(n0, n1, n2);
}

TypeDeterminer::PQState TypeDeterminer::create_pqstate(
											int n0, int n1, int n2) const {
	const double	p = genotype_probability(n0, n1, n2);
	return PQState(p, n0, n1, n2);
}

TypeDeterminer::PQState TypeDeterminer::create_pqstate(const State& s) const {
	const int	n0 = get<0>(s);
	const int	n1 = get<1>(s);
	const int	n2 = get<2>(s);
	return create_pqstate(n0, n1, n2);
}

TypeDeterminer::State TypeDeterminer::create_state(const PQState& s) const {
	const int	n0 = get<1>(s);
	const int	n1 = get<2>(s);
	const int	n2 = get<3>(s);
	return State(n0, n1, n2);
}

TypeDeterminer::State TypeDeterminer::reverse_state(const State& s) const {
	const int	n0 = get<0>(s);
	const int	n1 = get<1>(s);
	const int	n2 = get<2>(s);
	return State(n2, n1, n0);
}

vector<TypeDeterminer::PQState> TypeDeterminer::neighbor_states(
												const PQState& s0) const {
	vector<PQState>	neighbors;
	const int	n1 = get<1>(s0);
	const int	n2 = get<2>(s0);
	const int	n3 = get<3>(s0);
	if(n1 > 0) {
		if(n1-1 >= n3)
			neighbors.push_back(create_pqstate(n1-1, n2+1, n3));
		if(n1-1 >= n3+1)
			neighbors.push_back(create_pqstate(n1-1, n2, n3+1));
	}
	if(n2 > 0) {
		neighbors.push_back(create_pqstate(n1+1, n2-1, n3));
		if(n1 >= n3+1)
			neighbors.push_back(create_pqstate(n1, n2-1, n3+1));
	}
	if(n3 > 0) {
		neighbors.push_back(create_pqstate(n1+1, n2, n3-1));
		neighbors.push_back(create_pqstate(n1, n2+1, n3-1));
	}
	return neighbors;
}

// 0/1 x 0/1は難しい
void TypeDeterminer::make_memo11() {
	// 確率が大きい状態から並べて累積が1-αを超えるまで列挙する
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		const int	M = N - num_NA;
		const PQState	pqs0 = initialize_state(M);
		priority_queue<PQState>	pq;
		pq.push(pqs0);
		std::set<PQState>	visited;
		visited.insert(pqs0);
		double	total_p = 0.0;
		while(total_p < 1.0 - alpha) {
			const PQState	pqs = pq.top();
			pq.pop();
			if(get<1>(pqs) == get<3>(pqs)) {
				insert(create_state(pqs), 5);
				total_p += get<0>(pqs);
			}
			else {
				// n0 > n2なので、n0 < n2の分も考える
				const State	s1 = create_state(pqs);
				insert(s1, 5);
				insert(reverse_state(s1), 5);
				total_p += get<0>(pqs) * 2;
			}
			
			const vector<PQState>	neighbors = neighbor_states(pqs);
			for(auto p = neighbors.begin(); p != neighbors.end(); ++p) {
				const PQState	pqs1 = *p;
				if(visited.find(pqs1) == visited.end()) {
					pq.push(pqs1);
					visited.insert(pqs1);
				}
			}
		}
	}
}

void TypeDeterminer::make_memo12() {
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		for(int n0 = 0; n0 <= N / 10; ++n0) {
			const int	M = N - num_NA - n0;
			const int	n = binomial(M);
			for(int n1 = n; n1 <= M - n; ++n1)
				insert(n0, n1, M-n1, 6);
		}
	}
}

void TypeDeterminer::make_memo22() {
	// あり得ないGenotypeも1割だけあってもよいことにする
	for(int num_NA = 0; num_NA <= N / 10; ++num_NA) {
		for(int n0 = 0; n0 <= N / 10; ++n0) {
			for(int n1 = 0; n1 <= N / 10; ++n1)
				insert(n0, n1, N - n0 - n1 - num_NA, 10);
		}
	}
}

int TypeDeterminer::determine(const State& counter) const {
	const auto	p = memo.find(counter);
	if(p != memo.end() && p->second != OVERLAP)
		return p->second;
	else
		return -1;
}
