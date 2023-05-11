#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "../include/log.h"
#include "../include/BiasProbability.h"

using namespace std;

BiasProbability	*BiasProbability::bias_probability = NULL;


//////////////////// BiasProbability ////////////////////

int BiasProbability::compute_max_bias(int N, double float_cM) {
	const int	cM = (int)std::floor(float_cM);
	auto	p = ps.find(N);
	if(p == ps.end()) {
		ps[N].push_back(init(N));
		extend(ps[N], cM, N);
	}
	else if((double)ps[N].size() <= cM) {
		extend(ps[N], cM, N);
	}
	
	return compute_max_bias_core(ps[N][cM]);
}

BiasProbability::LineProbs BiasProbability::init(int N) {
	LineProbs	lps;
	for(int k = 0; k <= N; ++k) {
		auto	key = pair<int,int>(k, std::min(k, N-k));
		lps[key] = exp(log_C(N, k) + log(0.5) * N);
	}
	return lps;
}

void BiasProbability::extend(vector<LineProbs>& subps, int cM, int N) {
	for(int t = (int)subps.size(); t <= cM; ++t) {
		subps.push_back(step(subps[t-1], N));
	}
}

int BiasProbability::compute_max_bias_core(const LineProbs& lps) {
	map<int,double>	probs_max_bias;
	for(auto p = lps.begin(); p != lps.end(); ++p) {
		const int		max_bias = p->first.second;
		const double	prob = p->second;
		probs_max_bias[max_bias] += prob;
	}
	
	double	acc_prob = 0.0;
	for(int max_bias = 0; ; ++max_bias) {
		acc_prob += probs_max_bias[max_bias];
		if(acc_prob >= probability)
			return max_bias;
	}
}

BiasProbability::LineProbs BiasProbability::step(const LineProbs& lps, int N) {
	LineProbs	new_lps;
	for(auto p = lps.begin(); p != lps.end(); ++p) {
		const int	cur_bias = p->first.first;
		const int	max_bias = p->first.second;
		// numbers of 0->1, 1->0
		for(int i = 0; i <= N-cur_bias; ++i) {
			for(int j = 0; j <= cur_bias; ++j) {
				const int	new_bias = bias(cur_bias+i-j, N);
				const double	x = log_C(N-cur_bias, i) + log_C(cur_bias, j) +
										(i+j) * log(0.01) + (N-i-j) * log(0.99);
				const double	new_p = p->second * exp(x);
				const int		new_max_bias2 = std::min(max_bias, new_bias);
				new_lps[pair<int,int>(new_bias, new_max_bias2)] += new_p;
			}
		}
	}
	return new_lps;
}

// log(nCk)
double BiasProbability::log_C(int n, int k) {
	if(k == 0 || n == k)
		return 0.0;
	
	auto	p = memo_C.find(pair<int,int>(n, k));
	if(p != memo_C.end())
		return p->second;
	
	const double	res = Log::add(log_C(n-1, k-1), log_C(n-1, k));
	memo_C[pair<int,int>(n, k)] = res;
	return res;
}

int BiasProbability::bias(int i, int N) {
	return std::min(i, N-i);
}
