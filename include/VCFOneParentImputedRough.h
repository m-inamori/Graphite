#ifndef __VCFONEPARENTIMPUTEDROUGH
#define __VCFONEPARENTIMPUTEDROUGH

#include "VCFFamily.h"
#include "Map.h"

class Family;
class KnownFamily;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFOneParentImputedRough ////////////////////

class VCFOneParentImputedRough : public VCFBase, public VCFSmallBase,
								public VCFFamilyBase, public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const bool is_mat_imputed;
	const double	E[4][4];		// 排出確率
	const double	Epc[3][3][4];	// 両親と後代からの後代の排出確率
	
public:
	VCFOneParentImputedRough(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_imputed, const Map& map_, double w);
	~VCFOneParentImputedRough();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
	
private:
	///// common /////
	std::vector<int> trace_back(const std::vector<DP>& dps) const;
	
	///// parent /////
	int compute_non_phased_parent_gt(int h, int i) const;
	double emission_probability_for_parent(std::size_t i, int h, int op,
												const std::vector<int>& ocs,
												int phased_parent_gt) const;
	std::vector<DP> initialize_parent_dp() const;
	// hidden stateに対して、可能な前のhidden stateを集めておく
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	std::string get_phased_parent_gt(const VCFFamilyRecord *record) const {
		const size_t	phased_index = is_mat_imputed ? 0 : 1;
		return record->get_gt(phased_index);
	}
	std::string get_non_phased_parent_gt(const VCFFamilyRecord *record) const {
		const size_t	non_phased_index = is_mat_imputed ? 1 : 0;
		return record->get_gt(non_phased_index);
	}
	// genetic distance between two records
	double dist(const VCFRecord *r1, const VCFRecord *r2) const;
	void update_dp_for_parent(std::size_t i, std::vector<DP>& dp,
					const std::vector<std::vector<int>>& prev_h_table) const;
	void set_non_phased_parent_gt(int gt, VCFFamilyRecord *record);
	void update_parent_genotypes(const std::vector<int>& hs);
	void impute_parent();
	
	///// progenies /////
	double emission_probability_for_progeny(int h, int oc,
											int mat_gt, int pat_gt) const;
	std::vector<DP> initialize_progeny_dp(std::size_t ic) const;
	void update_dp_for_progeny(std::size_t i, std::size_t ic,
												std::vector<DP>& dp) const;
	void update_progeny_genotypes(const std::vector<int>& hs, int ic);
	void impute_progeny(std::size_t ic);
	
private:
	static int gt_by_haplotypes(int hc1, int hc2, int mat_gt, int pat_gt);
};
#endif
