#ifndef __VCFLARGEBOTHREF
#define __VCFLARGEBOTHREF

#include "VCF.h"
#include "VCFImputable.h"
#include "ProgenyRefImputer.h"

class Map;


//////////////////// VCFLargeBothRef ////////////////////

class VCFLargeBothRef : public VCFImputable {
	std::vector<VCFFamilyRecord *>	records;
	ProgenyRefImputer	imputer;
	
public:
	VCFLargeBothRef(const STRVEC& samples,
					const std::vector<VCFFamilyRecord *>& rs,
					const Map& map_, double w, const VCFSmall *vcf):
										VCFImputable(samples, vcf),
										records(rs),
										imputer(records, map_, w) { }
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	std::size_t amount() const override { return 1; }
	void impute() override;
};

#endif
