#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "../include/log.h"
#include "../include/Baum_Welch_with_fixed_Ts.h"

using namespace std;


BaumWelch::BaumWelch(const string& seq_,
						const vector<char>& hidden_states_,
						const vector<char>& states_,
						const vector<Matrix>& Ts_) :
								hidden_states(hidden_states_), states(states_),
								seq(seq_),
								Ts(modify_log(Ts_)),
								pi(initialize_pi()),
//								T(initialize_transition_matrix()),
								E(initialize_emission_matrix()) {
	if(!converge_pi_E())
		converge_pi();
}

vector<BaumWelch::Matrix> BaumWelch::modify_log(const vector<Matrix>& Ms_) {
	vector<Matrix>	Ms;
	for(auto p = Ms_.begin(); p != Ms_.end(); ++p) {
		Matrix	M;
		for(auto q = p->begin(); q != p->end(); ++q)
			M[q->first] = Log::modified_log(q->second);
		Ms.push_back(M);
	}
	return Ms;
}

map<char,double> BaumWelch::initialize_pi() {
	const double	prob = -std::log((double)hidden_states.size());
	map<char,double>	probs;
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p)
		probs[*p] = prob;
	return probs;
}

BaumWelch::Matrix BaumWelch::initialize_transition_matrix() {
	Matrix	T;
	const size_t	N = hidden_states.size();
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		for(auto q = hidden_states.begin(); q != hidden_states.end(); ++q) {
			set(T, *p, *q, std::log(*p == *q ? 0.9 : 0.1 / (N - 1)));
		}
	}
	return T;
}

BaumWelch::Matrix BaumWelch::initialize_emission_matrix() {
	// probability of hidden state -> state
	Matrix	E;
	const size_t	N = states.size();
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		for(auto q = states.begin(); q != states.end(); ++q) {
			set(E, *p, *q, std::log(*p == *q ? 0.9 : 0.1 / (N - 1)));
		}
	}
	return E;
}

BaumWelch::Table BaumWelch::compute_alpha(const Matrix& A) const {
	const size_t	L = seq.size();
	Table	a(L);	// probability of s
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		auto	q = pi.find(*p);
		assert(q != pi.end());
		a[0][*p] = q->second + get(A, seq.c_str()[0], *p);
	}
	for(size_t t = 1U; t < seq.size(); ++t) {
		const Matrix&	T = Ts[t-1];
		const char		s = seq.c_str()[t];
		for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
			const char	h = *p;
			vector<double>	v;
			for(auto q = hidden_states.begin(); q != hidden_states.end(); ++q) {
				const char	h_prev = *q;
				v.push_back(a[t-1][h_prev] + get(T, h_prev, h) + get(A, s, h));
			}
			a[t][h] = Log::sum(v);
		}
	}
	return a;
}

BaumWelch::Table BaumWelch::compute_beta(const Matrix& A) const {
	const size_t	L = seq.size();
	Table	b(L);
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p)
		b[L-1][*p] = 0.0;	// log(1.0)
	for(int t = (int)L - 2; t >= 0; --t) {
		const Matrix&	T = Ts[t];
		const char		s = seq.c_str()[t+1];
		for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
			const char	h = *p;
			vector<double>	v;
			for(auto q = hidden_states.begin(); q != hidden_states.end(); ++q) {
				const char	h_next = *q;
				v.push_back(get(T, h, h_next) + get(A, s, h_next) +
														b[t+1][h_next]);
			}
			b[t][h] = Log::sum(v);
		}
	}
	return b;
}

BaumWelch::Table BaumWelch::compute_gamma(
						const Table& a, const Table& b) const {
	Table	gamma1;
	for(size_t t = 0; t < a.size(); ++t) {
		auto&	a1 = a[t];
		auto&	b1 = b[t];
		map<char,double>	dic;
		for(auto r = hidden_states.begin(); r != hidden_states.end(); ++r) {
			const char	h = *r;
			auto	pa = a1.find(h);
			auto	pb = b1.find(h);
			assert(pa != a1.end() && pb != b1.end());
			dic[h] = pa->second + pb->second;
//			dic[h] = a1[h] + b1[h];
		}
		gamma1.push_back(dic);
	}
	
	double	prob_seq;	// probability that seq is observed
	{
		vector<double>	v;
		for(auto p = gamma1[0].begin(); p != gamma1[0].end(); ++p)
			v.push_back(p->second);
		prob_seq = Log::sum(v);
	}
	
	Table	gamma;
	for(auto p = gamma1.begin(); p != gamma1.end(); ++p) {
		map<char,double>	dic;
		const map<char,double>&	g = *p;
		for(auto q = g.begin(); q != g.end(); ++q)
			dic[q->first] = q->second - prob_seq;
		gamma.push_back(dic);
	}
	
	return gamma;
}

BaumWelch::Table BaumWelch::compute_parameters(const Matrix& A) {
	// probability of s
	const Table	a = compute_alpha(A);
	const Table	b = compute_beta(A);
	return compute_gamma(a, b);
}

void BaumWelch::update_probabilities(map<char,double>& pi1, Matrix& E1) {
	const Matrix	A = reverse_probs(E);
	Table	gamma = compute_parameters(A);
	pi1 = gamma[0];
	
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		const char	hj = *p;
		vector<double>	v;
		for(auto p = gamma.begin(); p != gamma.end(); ++p) {
			auto	g = *p;
			v.push_back(g[hj]);
		}
		const double	cum = Log::sum(v);
		for(auto q = states.begin(); q != states.end(); ++q) {
			const char	s = *q;
			vector<double>	w;
			for(size_t t = 0U; t < seq.size(); ++t) {
				const char	s1 = seq.c_str()[t];
				auto	g = gamma[t];
				if(s1 == s)
					w.push_back(g[hj]);
			}
			if(w.empty())
				continue;
			set(E1, hj, s, Log::sum(w) - cum);
		}
	}
}

bool BaumWelch::is_valid_Emission_Matrix(const Matrix& E) const {
	for(auto p = E.begin(); p != E.end(); ++p) {
		const auto	state = p->first;
		const double	v = p->second;
		if(state.first == state.second && v <= std::log(0.5))
			return false;
	}
	return true;
}

bool BaumWelch::is_near_parameters(const Matrix& E1,
								   const map<char,double>& pi1,
								   const Matrix& E2,
								   const map<char,double>& pi2) const {
	vector<double>	v1;
	for(auto p = E1.begin(); p != E1.end(); ++p) {
		auto	q = E2.find(p->first);
		if(q != E2.end())
			v1.push_back(Log::diff(p->second, q->second));
		else
			v1.push_back(Log::LOGZERO);
	}
	const double	d2 = Log::sum(v1);
	
	vector<double>	v2;
	for(auto p = pi1.begin(); p != pi1.end(); ++p) {
		auto	q = pi2.find(p->first);
		assert(q != pi2.end());
		v2.push_back(Log::diff(p->second, q->second));
	}
	const double	d3 = Log::sum(v2);
	
	return Log::add(d2, d3) < std::log(1e-7);
}

bool BaumWelch::converge_pi_E() {
	for(int i = 0; i < 100; ++i) {
		try {
			map<char,double>	pi1;
			Matrix				E1;
			update_probabilities(pi1, E1);
			if(is_near_parameters(E, pi, E1, pi1))
				break;
			pi = pi1;
			E =  E1;
		}
		catch(runtime_error& e) {
			break;
		}
	}
	
	return is_valid_Emission_Matrix(E);
}

BaumWelch::Matrix BaumWelch::reverse_probs(const Matrix& M) const {
	map<char,double>	sum_observed;
	for(auto p = states.begin(); p != states.end(); ++p) {
		const char	s0 = *p;
		vector<double>	v;
		for(auto q = M.begin(); q != M.end(); ++q) {
			const char		s = q->first.second;
			const double	p = q->second;
			if(s == s0)
				v.push_back(p);
		}
		sum_observed[s0] = Log::sum(v);
	}
	
	Matrix	rev_M;
	for(auto p = M.begin(); p != M.end(); ++p) {
		const char	h = p->first.first;
		const char	s = p->first.second;
		auto	q = M.find(p->first);
		assert(q != M.end());
		rev_M[pair<char,char>(s, h)] = q->second - sum_observed[s];
	}
	return rev_M;
}

map<char,double> BaumWelch::update_pi(map<char,double>& pi) {
	const auto	A = reverse_probs(E);
	const auto	gamma = compute_parameters(A);
	const auto	pi1 = gamma.front();
	return pi1;
}

bool BaumWelch::converge_pi() {
	pi = initialize_pi();
	E = initialize_emission_matrix();
	for(int i = 0; i < 100; ++i) {
		try {
			auto	pi1 = update_pi(pi);
			vector<double>	v;
			for(auto p = pi.begin(); p != pi.end(); ++p)
				v.push_back(Log::diff(p->second, pi1[p->first]));
			if(Log::sum(v) < std::log(1e-7))
				break;
			pi = pi1;
		}
		catch(runtime_error& e) {
			break;
		}
	}
	return true;
}

string BaumWelch::Viterbi(const Matrix& A) const {
	Table	table(seq.size());
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		const char	h = *p;
		auto	q = pi.find(h);
		assert(q != pi.end());
		table[0][h] = q->second + get(A, seq[0], h);
	}
	
	for(size_t k = 1U; k < seq.size(); ++k) {
		const char	s = seq[k];
		const auto&	T = Ts[k-1];
		for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
			const char	h = *p;
			double	max_value = Log::LOGZERO;
			for(auto q = hidden_states.begin(); q != hidden_states.end(); ++q) {
				const char	h0 = *q;
				const double	v = table[k-1][h0] + get(T, h0,h) + get(A, s,h);
				max_value = std::max(max_value, v);
			}
			table[k][h] = max_value;
		}
	}
	
	// backtrack
	double	prob = Log::LOGZERO;
	char	h = '0';	// temporary
	string	new_hidden_seq = "";
	for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
		const char	h1 = *p;
		double		v = table.back()[h1];
		if(v > prob) {
			prob = v;
			h = h1;
		}
	}
	new_hidden_seq += h;
	
	for(size_t k = seq.size() - 1; k > 0U; --k) {
		const char	s = seq.c_str()[k];
		const Matrix&	T = Ts[k-1];
		char	h1 = '0';	// temporary
		double	prob1 = Log::LOGZERO;
		for(auto p = hidden_states.begin(); p != hidden_states.end(); ++p) {
			const char	h0 = *p;
			double	v = table[k-1][h0] + get(T, h0, h) + get(A, s, h);
			if(v > prob1) {
				prob1 = v;
				h1 = h0;
			}
		}
		new_hidden_seq += h1;
		h = h1;
	}
	
	std::reverse(new_hidden_seq.begin(), new_hidden_seq.end());
	return new_hidden_seq;
}

BaumWelch::Matrix BaumWelch::reverse_emission_matrix() const {
	return this->reverse_probs(E);
}

double BaumWelch::get(const Matrix& M, char s1, char s2) {
	auto	p = M.find(pair<char,char>(s1, s2));
	assert(p != M.end());
	return p->second;
}

void BaumWelch::set(Matrix& M, char s1, char s2, double v) {
	M.insert(make_pair(pair<char,char>(s1, s2), v));
}

string BaumWelch::impute(const string& seq,
							const vector<char> hidden_states,
							const vector<char> states,
							const vector<Matrix>& Ts) {
	BaumWelch	BW(seq, hidden_states, states, Ts);
	const auto	A = BW.reverse_emission_matrix();
	return BW.Viterbi(A);
}
