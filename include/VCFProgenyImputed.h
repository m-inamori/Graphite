#ifndef __VCFPROGENYIMPUTED
#define __VCFPROGENYIMPUTED

#include "VCFFamily.h"
#include "ParentImputerByProgeny.h"

class VCFSmall;


//////////////////// VCFProgenyImputed ////////////////////

class VCFProgenyImputed : public VCFFamilyBase {
	const std::vector<VCFFamilyRecord *>	records;
	ParentImputerByProgeny	imputer;
	
public:
	VCFProgenyImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool is_mat_known, const Map& map_,
						double w, const VCFSmall *vcf);
	VCFProgenyImputed(const VCFProgenyImputed&) = delete;
	VCFProgenyImputed& operator=(const VCFProgenyImputed&) = delete;
	~VCFProgenyImputed();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute();
};
#endif
