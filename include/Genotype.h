#ifndef __GENOTYPE
#define __GENOTYPE

#include <string>
#include <vector>
#include <utility>


//////////////////// Genotype ////////////////////

class Genotype {
	const char	gt1;
	const char	gt2;
	const bool	phasing;
	
public:
	explicit Genotype(const std::string& s);
	
	std::pair<char,char> gts() const;
	bool includes(char gt) const;
	bool conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
													bool considers_phasing);
	
	static bool is_homo(const std::string& gt) {
		return gt.length() >= 3 &&
				(gt.c_str()[0] == '0' or gt.c_str()[0] == '1') &&
				gt[0] == gt[2];
	}
	static int get_int_gt(const std::string& s);
	static std::string int_to_gt(int n);
	static bool is_NA(const std::string& gt) {
		return gt.c_str()[0] == '.' || gt.c_str()[2] == '.';
	}
	static bool is_valid(const std::string& gt, int mat_gt, int pat_gt);
	static std::string possible_gts(int gt);
	static int sum_gt(const std::string& gt);
	static bool is_all_NA(const std::vector<std::string>& GTs);
	static std::string int_to_phased_gt(int gt_int);
	
	static int phased_gt_to_int(const std::string& gt) {
		const int	gt1_ = gt.c_str()[0] == '0' ? 0 : 1;
		const int	gt2_ = gt.c_str()[2] == '0' ? 0 : 1;
		return gt1_ | (gt2_ << 1);
	}
	// 上のget_int_gtとN/Aのときが違う
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
