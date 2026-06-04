#ifndef __VCFPROGENYIMPUTED
#define __VCFPROGENYIMPUTED

#include "VCFImputable.h"
#include "ParentImputerByProgeny.h"

class VCFSmall;


//////////////////// VCFProgenyImputed ////////////////////

class VCFProgenyImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const bool	should_impute_mat;
	ParentImputerByProgeny	imputer;
	
public:
	VCFProgenyImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						bool should_impute_mat, const Map& map_,
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
	
	///// virtual methods for VCFFamilyBase /////
	std::size_t amount() const override { return 1; }
	STRVEC imputed_samples() const override;
	void impute() override;
};
#endif
