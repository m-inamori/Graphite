#ifndef __VCFBOTHPARENTIMPUTED
#define __VCFBOTHPARENTIMPUTED

#include "VCFImputable.h"
#include "ProgenyImputer.h"

class Map;


//////////////////// VCFBothParentImputed ////////////////////

class VCFBothParentImputed : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	ProgenyImputer	prog_imputer;
	
public:
	VCFBothParentImputed(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& records,
							const Map& map_, double w, const VCFSmall *vcf);
	VCFBothParentImputed(const VCFBothParentImputed&) = delete;
	VCFBothParentImputed& operator=(const VCFBothParentImputed&) = delete;
	~VCFBothParentImputed();
	
	///// virtual methods for VCFGenoBase /////
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	std::size_t size() const override { return records.size(); }
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	std::size_t amount() const override { return 1; }
	void impute() override;
};
#endif
