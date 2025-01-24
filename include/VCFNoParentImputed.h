#ifndef __VCFNOPARENTIMPUTED
#define __VCFNOPARENTIMPUTED

#include "VCFFamily.h"
#include "Map.h"

class Family;
class KnownFamily;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFNoParentImputed ////////////////////

class VCFNoParentImputed : public VCFBase, public VCFSmallBase,
								public VCFFamilyBase, public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	const double	MIN_PROB = -1e300;
	
private:
	std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const double	E[4][4];			// exhaust probabilities
	const double	Epc[3][4][4];		// progeny's exhaust probabilities
	const std::vector<double>	Cp;		// crossover values for parent
	
public:
	VCFNoParentImputed(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w);
	~VCFNoParentImputed();
	
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
	void clear_records() { records.clear(); }
	
private:
	///// common /////
	std::vector<int> trace_back(const std::vector<DP>& dps) const;
	
	///// parent /////
	int compute_non_phased_parent_gt(int h, int i) const;
	double emission_probability(std::size_t i, int h,
										int mat_gt, int pat_gt) const;
	std::vector<DP> initialize_dp() const;
	// collect the possible previous hidden states for each hidden state
	std::vector<std::vector<int>>
				collect_possible_previous_hidden_states() const;
	// genetic distance between two records
	double dist(const VCFRecord *r1, const VCFRecord *r2) const;
	std::vector<double> calc_Cp(const std::vector<VCFFamilyRecord *>& rs) const;
	void update_dp(std::size_t i, std::vector<DP>& dp,
					const std::vector<std::vector<int>>& prev_h_table) const;
	void set_phased_gt(int gt, VCFFamilyRecord *record);
	void update_genotypes(const std::vector<int>& hs);
	void impute_parent();
};
#endif
