#ifndef __GENOTYPE
#define __GENOTYPE

#include <string>
#include <vector>
#include <utility>


//////////////////// Genotype ////////////////////

class Genotype {
public:
	static const int	UN_00 = 0;	// 0/0
	static const int	UN_01 = 1;	// 0/1
	static const int	UN_11 = 2;	// 1/1
	static const int	NA = 3;		// ./.
	static const int	PH_00 = 4;	// 0|0
	static const int	PH_10 = 5;	// 1|0
	static const int	PH_01 = 6;	// 0|1
	static const int	PH_11 = 7;	// 1|1
	
	const char	gt1;
	const char	gt2;
	const bool	phasing;
	
public:
	explicit Genotype(const std::string& s);
	
	std::pair<char,char> gts() const;
	bool includes(char gt) const;
	bool conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
													bool considers_phasing);
	
	static bool is_ref_homo(int gt) {
		return gt == Genotype::UN_00 || gt == Genotype::PH_00;
	}
	static bool is_alt_homo(int gt) {
		return gt == Genotype::UN_11 || gt == Genotype::PH_11;
	}
	static bool is_homo(int gt) {
		if(gt < 6)
			return (gt & 1) == 0;
		else
			return gt == 7;
	}
	static bool is_hetero(int gt) {
		return gt == Genotype::UN_01 || gt == Genotype::PH_10 ||
										gt == Genotype::PH_01;
	}
	
	static int unphased(int gt) {
		if(gt < 4)			return gt;
		else if(gt == 4)	return 0;
		else if(gt < 7)		return 1;
		else				return 2;
	}
	
	static std::string int_to_gt(int n);
	static bool is_valid(int gt, int mat_gt, int pat_gt);
	static std::vector<int> possible_gts(int gt);
	
	static int phased_gt_to_int(const std::string& gt) {
		const int	gt1_ = gt.c_str()[0] == '0' ? 0 : 1;
		const int	gt2_ = gt.c_str()[2] == '0' ? 0 : 1;
		return gt1_ | (gt2_ << 1);
	}
	
	// N/A: 3
	// non-phased: 0-2
	// phased: 4-7
	static int all_gt_to_int(const std::string& gt) {
		if(gt.length() <= 2 || gt.c_str()[0] == '.' || gt.c_str()[2] == '.')
			return Genotype::NA;
		else if(gt.c_str()[1] == '/')
			return Genotype::gt_to_int(gt);
		else
			return Genotype::phased_gt_to_int(gt) | 4;
	}
	
	// 0-7 -> 0-3
	static int all_int_gt_to_int_gt(int gt_int) {
		if(gt_int < 5)
			return gt_int & 3;
		else
			// 6(0|1) -> 1
			// 7(1|1) -> 2
			return gt_int - 5;
	}
	
	// two alleles -> phased genotype
	static int from_alleles(int a1, int a2) {
		return a1 | (a2 << 1) | 4;
	}
	
	static bool is_NA(int gt) {
		return gt == Genotype::NA;
	}
	
	static bool is_00(int gt) {
		return (gt & 3) == 0;
	}
	
	static bool is_01(int gt) {
		return gt == Genotype::UN_01 ||
			   gt == Genotype::PH_01 || gt == Genotype::PH_10;
	}
	
	static bool is_11(int gt) {
		return gt == Genotype::UN_11 or gt == Genotype::PH_11;
	}
	
	static bool is_phased(int gt) {
		return gt >= 4;
	}
	
	static std::string int_to_phased_gt(int gt) {
		switch(gt & 3) {
			case  0: return "0|0";
			case  1: return "1|0";
			case  2: return "0|1";
			default: return "1|1";
		}
	}
	
	static std::string int_to_all_gt(int gt_int) {
		if(gt_int < 4)
			return Genotype::int_to_gt(gt_int);
		else
			return Genotype::int_to_phased_gt(gt_int & 3);
	}
	
	// 0|1 <-> 1|0 provide phased
	static int inverse(int gt) {
		switch(gt) {
			case Genotype::PH_01: return Genotype::PH_10;
			case Genotype::PH_10: return Genotype::PH_01;
			default:              return gt;
		}
	}
	
	// allele of j(0 or 1) side of phased Genotype
	static int get_allele(int gt, int j) {
		return (gt >> j) & 1;
	}
	
	static int inverse_allele(int a, bool inv) {
		if(inv)
			return a == 0 ? 1 : 0;
		else
			return a;
	}
	
	static int gt_by_haplotypes(int hc1, int hc2, int mat_gt, int pat_gt) {
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1);
	}
	
	static int gt_to_int(const std::string& gt) {
		for(int i = 0; i < 3; ++i) {
			switch(gt.c_str()[i]) {
				case '.':  return 3;
				case '\0': return 3;
			}
		}
		const int	gt1_ = gt.c_str()[0] == '0' ? 0 : 1;
		const int	gt2_ = gt.c_str()[2] == '0' ? 0 : 1;
		return gt1_ + gt2_;
	}
	static std::size_t find_key_position(const std::string& info,
											const std::string& key);
};
#endif
