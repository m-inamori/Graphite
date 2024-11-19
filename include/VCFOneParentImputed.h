#ifndef __VCFONEPARENTIMPUTED
#define __VCFONEPARENTIMPUTED

#include "VCFFamily.h"
#include "Map.h"

class Family;
class KnownFamily;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFOneParentImputed ////////////////////

class VCFOneParentImputed : public VCFBase, public VCFSmallBase,
							public VCFFamilyBase, public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>
	
private:
	const std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	const bool is_mat_imputed;
	const float	E[4][4];		// 排出確率
	
	VCFOneParentImputed(const std::vector<STRVEC>& header,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>& ref_haps,
						bool is_mat_imputed, const Map& map_) :
										VCFBase(header), VCFSmallBase(),
										VCFFamilyBase(), VCFMeasurable(map_)
	~VCFOneParentImputed() { }
	
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
	VCFFamilyRecord *get_family_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
