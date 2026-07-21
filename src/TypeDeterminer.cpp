#include <set>
#include <queue>
#include <algorithm>
#include <cmath>
#include "../include/TypeDeterminer.h"
#include "../include/MathCommon.h"

using namespace std;


//////////////////// TypeDeterminer ////////////////////

vector<double> TypeDeterminer::make_memo_log_facs() {
	const size_t	L = N + b + 3;
	vector<double>	v(L, 0.0);
	for(size_t n = 1; n < L; ++n) {
		v[n] = v[n-1] + log(static_cast<double>(n));
	}
	return v;
}

// in the case which all genotypes are the same without errors
double TypeDeterminer::log_likelihood_one(size_t n, size_t m) const {
	return log_beta(m + a, n + b) - m*log(3.0) - log_beta(a, b);
}

// log likelihood for 0/0 x 0/0
double TypeDeterminer::log_likelihood00(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 0/0 x 0/0 -> 0/0: 1-p otherwise: p/3
	const int	n = gt_freq[0];								// 0/0
	const int	m = gt_freq[1] + gt_freq[2] + gt_freq[3];	// 0/1 + 1/1 + ./.
	if(mat_gt == 0 && pat_gt == 0) {
		// likelihood of both parents: (1-p)^2
		return log_likelihood_one(n+2, m);
	}
	else if(mat_gt == 0 || pat_gt == 0) {
		// only one parent matches
		// likelihood of both parents: (1-p)p
		return log_likelihood_one(n+1, m+1);
	}
	else {
		// neither parent matches
		// likelihood of both parents: p^2
		return log_likelihood_one(n, m+2);
	}
}

// log likelihood for 0/0 x 0/1
double TypeDeterminer::log_likelihood01(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 0/0 x 0/1 ->	-> 0/0 or 0/1: 1/2-p/3, otherwise: p/3
	const int	n = gt_freq[0] + gt_freq[1];	// 0/0 or 0/1
	const int	m = gt_freq[2] + gt_freq[3];	// 1/1 or ./.
	if((mat_gt == 0 && pat_gt == 1) || (mat_gt == 1 && pat_gt == 0)) {
		// likelihood of both parents: (1-p)^2
		auto f = [n, m](double p) {
					return std::pow(1.0 - p, 2)
						 * std::pow(0.5 - p / 3.0, n) 
						 * std::pow(p / 3.0, m) 
						 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else if(mat_gt == 0 || mat_gt == 1 || pat_gt == 0 || pat_gt == 1) {
		// only one parent matches
		// likelihood of both parents: (1-p)p/3
		auto f = [n, m](double p) {
					return (1.0 - p) * p / 3.0 
						 * std::pow(0.5 - p / 3.0, n) 
						 * std::pow(p / 3.0, m) 
						 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else {
		// both parents not matche
		// likelihood of both parents: (p/3)^2
		auto f = [n, m](double p) {
					return std::pow(p / 3.0, 2)
						 * std::pow(0.5 - p / 3.0, n) 
						 * std::pow(p / 3.0, m) 
						 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
}

// log likelihood for 0/1 x 0/1
double TypeDeterminer::log_likelihood11(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 0/1 x 0/1 -> 0/0: 1/4, 0/1: 1/2-p/3, 1/1: 1/4, otherwise: p/3
	const int	n = gt_freq[0] + gt_freq[2];
	const int	m = gt_freq[1];
	const int	q = gt_freq[3];
	
	if(mat_gt == 1 && pat_gt == 1) {
		// both parents match
		// likelihood for both parents is (1-p)^2
		auto f = [n, m, q](double p) {
			return std::pow(1.0 - p, 2) 
				 * std::pow(1.0 / 4.0, n) 
				 * std::pow(1.0 / 2.0 - p / 3.0, m) 
				 * std::pow(p / 3.0, q) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else if(mat_gt == 1 || pat_gt == 1) {
		// only one parent matches
		// likelihood for both parents is (1-p)p/3
		auto f = [n, m, q](double p) {
			return (1.0 - p) * p / 3.0 
				 * std::pow(1.0 / 4.0, n) 
				 * std::pow(1.0 / 2.0 - p / 3.0, m) 
				 * std::pow(p / 3.0, q) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else {
		// neither parent matches
		// likelihood for both parents is (p/3)^2
		auto f = [n, m, q](double p) {
			return std::pow(p / 3.0, 2) 
				 * std::pow(1.0 / 4.0, n) 
				 * std::pow(1.0 / 2.0 - p / 3.0, m) 
				 * std::pow(p / 3.0, q) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
    }
}

// log likelihood for 0/0 x 1/1
double TypeDeterminer::log_likelihood02(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 0/0 x 1/1 ->	-> 0/1: 1-p, otherwise: p/3
	const int	n = gt_freq[1];
	const int	m = gt_freq[0] + gt_freq[2] + gt_freq[3];
	
	if((mat_gt == 0 && pat_gt == 2) || (mat_gt == 2 && pat_gt == 0)) {
		// both parents match
		// likelihood for both parents is (1-p)^2
		return log_likelihood_one(n + 2, m);
	}
	else if(mat_gt == 0 || mat_gt == 2 || pat_gt == 0 || pat_gt == 2) {
		// only one parent matches
		// likelihood for both parents is (1-p)p/3
		return log_likelihood_one(n + 1, m + 1);
	}
	else {
		// neither parent matches
		// likelihood for both parents is (p/3)^2
		return log_likelihood_one(n, m + 2);
	}
}

// log likelihood for 0/1 x 1/1
double TypeDeterminer::log_likelihood12(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 0/1 x 1/1 -> 0/1 or 1/1: 1/2-p/3, otherwise: p/3
	const int	n = gt_freq[1] + gt_freq[2];	// 0/1 or 1/1
	const int	m = gt_freq[0] + gt_freq[3];	// 0/0 or ./.
	
	if((mat_gt == 1 && pat_gt == 2) || (mat_gt == 2 && pat_gt == 1)) {
		// both parents match
		// likelihood for both parents is (1-p)^2
		auto f = [n, m](double p) {
			return std::pow(1.0 - p, 2) 
				 * std::pow(1.0 / 2.0 - p / 3.0, n) 
				 * std::pow(p / 3.0, m) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else if(mat_gt == 1 || mat_gt == 2 || pat_gt == 1 || pat_gt == 2) {
		// only one parent matches
		// likelihood for both parents is (1-p)p/3
		auto f = [n, m](double p) {
			return (1.0 - p) * p / 3.0 
				 * std::pow(1.0 / 2.0 - p / 3.0, n) 
				 * std::pow(p / 3.0, m) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
	else {
		// neither parent matches
		// likelihood for both parents is p^2
		auto f = [n, m](double p) {
			return std::pow(p / 3.0, 2) 
				 * std::pow(1.0 / 2.0 - p / 3.0, n) 
				 * std::pow(p / 3.0, m) 
				 * std::pow(1.0 - p, 8);
		};
		return log(MathCommon::Simpson(f, 0.0, 1.0, 10)) - log_beta(a, b);
	}
}

// log likelihood for 1/1 x 1/1
double TypeDeterminer::log_likelihood22(int mat_gt, int pat_gt,
										const array<int, 4>& gt_freq) const {
	// 1/1 x 1/1 -> 1/1: 1-p otherwise: p/3
	const int	n = gt_freq[2];								// 1/1
	const int	m = gt_freq[0] + gt_freq[1] + gt_freq[3];	// 0/0 + 0/1 + ./.
	
	if(mat_gt == 2 && pat_gt == 2) {
		// both parents match
		// likelihood for both parents is (1-p)^2
		return log_likelihood_one(n + 2, m);
	}
	else if(mat_gt == 2 || pat_gt == 2) {
		// only one parent matches
		// likelihood for both parents is (1-p)p/3
		return log_likelihood_one(n + 1, m + 1);
	}
	else {
		// neither parent matches
		// likelihood for both parents is (p/3)^2
		return log_likelihood_one(n, m + 2);
	}
}

vector<pair<ParentComb, double>> TypeDeterminer::determine(
										int mat_gt, int pat_gt,
										const array<int, 4>& counter) const {
	const double	l00 = log_likelihood00(mat_gt, pat_gt, counter);
	const double	l01 = log_likelihood01(mat_gt, pat_gt, counter);
	const double	l11 = log_likelihood11(mat_gt, pat_gt, counter);
	const double	l02 = log_likelihood02(mat_gt, pat_gt, counter);
	const double	l12 = log_likelihood12(mat_gt, pat_gt, counter);
	const double	l22 = log_likelihood22(mat_gt, pat_gt, counter);
	
	const double	ls[6] = { l00, l01, l11, l02, l12, l22 };
	
	int		max_index = 0;
	double	max_l = ls[0];
	for (int i = 1; i < 6; ++i) {
		if (ls[i] > max_l) {
			max_l = ls[i];
			max_index = i;
		}
	}
	
	std::vector<std::pair<ParentComb, double>>	pairs;
	pairs.push_back(std::make_pair(static_cast<ParentComb>(max_index), max_l));
	
	const double	threshold = max_l - std::log(30.0);
	for (int i = 0; i < 6; ++i) {
		if (i != max_index && ls[i] > threshold) {
			pairs.push_back(std::make_pair(static_cast<ParentComb>(i), ls[i]));
		}
	}
	
	return pairs;
}
