#ifndef __VCFPARENTSKNOWNPROGENYIMPUTED
#define __VCFPARENTSKNOWNPROGENYIMPUTED

#include "VCFImputable.h"
#include "ParentsImputerByProgeny.h"


//////////////////// VCFParentsKnownProgenyImputed ////////////////////

class VCFParentsKnownProgenyImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const std::size_t	NH;
	ParentsImputerByProgeny	imputer;
	
public:
	VCFParentsKnownProgenyImputed(const STRVEC& s,
						const std::vector<VCFFamilyRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps_mat,
						const std::vector<std::vector<int>>& ref_haps_pat,
						const Map& map_, double w, const VCFSmall *vcf);
	VCFParentsKnownProgenyImputed(
					const VCFParentsKnownProgenyImputed&) = delete;
	VCFParentsKnownProgenyImputed&
					operator=(const VCFParentsKnownProgenyImputed&) = delete;
	~VCFParentsKnownProgenyImputed();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	void impute() override;
	STRVEC imputed_samples() const override {
		return { samples[0], samples[1] };
	}
	std::size_t amount() const override;
};
#endif
