#ifndef __VCFPROGENYIMPUTED
#define __VCFPROGENYIMPUTED

#include "VCF.h"
#include "Map.h"

class ParentImputerByProgeny;


//////////////////// VCFProgenyImputed ////////////////////

class VCFProgenyImputed : public VCFBase, public VCFSmallBase {
	const std::vector<VCFRecord *>	records;
	ParentImputerByProgeny	*imputer;
	
public:
	VCFProgenyImputed(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_known, const Map& map_, double w);
	VCFProgenyImputed(const VCFProgenyImputed&) = delete;
	VCFProgenyImputed& operator=(const VCFProgenyImputed&) = delete;
	~VCFProgenyImputed();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
