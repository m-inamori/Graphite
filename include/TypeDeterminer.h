#ifndef __TYPEDETERMINER
#define __TYPEDETERMINER

#include <vector>
#include <map>
#include <tuple>
#include <array>


//////////////////// ParentComb ////////////////////

enum class ParentComb {
	P00x00 = 0, P00x01, P01x01, P00x11, P01x11, P11x11, PNA
};


//////////////////// TypeDeterminer ////////////////////

class TypeDeterminer {
	const std::size_t	N;
	const double	a;
	const double	b;
	const std::vector<double>	log_facs;
	
public:
	TypeDeterminer(size_t n) : N(n), a(1), b(9),
								log_facs(make_memo_log_facs()) { }
	
	std::vector<std::pair<ParentComb, double>> determine(
										int mat_gt, int pat_gt,
										const std::array<int, 4>& counter) const;
	
private:
	std::vector<double> make_memo_log_facs();
	
	double log_beta(std::size_t n, std::size_t m) const {
		return log_facs[n-1] + log_facs[m-1] - log_facs[n+m-1];
	}
	
	// Cases where all result in the same genotype in the absence of errors
	double log_likelihood_one(std::size_t n, std::size_t m) const;
	// log likelihood for 0/0 x 0/0
	double log_likelihood00(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	// log likelihood for 0/0 x 0/1
	double log_likelihood01(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	// log likelihood for 0/1 x 0/1
	double log_likelihood11(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	// log likelihood for 0/0 x 1/1
	double log_likelihood02(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	// log likelihood for 0/1 x 1/1
	double log_likelihood12(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	// log likelihood for 1/1 x 1/1
	double log_likelihood22(int mat_gt, int pat_gt,
							const std::array<int, 4>& gt_freq) const;
	
public:
	// Since an enum does not act as a namespace,
	// I am placing the functions for ParentComb here.
	static bool is_same_parent_genotype(ParentComb c) {
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
	
	static std::pair<int,int> int_gt_pair(ParentComb comb) {
		const int	p = static_cast<int>(comb);
		if(p == 0)
			return std::pair<int, int>(0, 0);
		else if(p < 3)
			return std::pair<int, int>(p - 1, 1);
		else
			return std::pair<int, int>(p - 3, 2);
	}
};
#endif
