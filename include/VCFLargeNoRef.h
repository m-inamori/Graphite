#ifndef __VCFLARGENOREF
#define __VCFLARGENOREF

#include "VCF.h"
#include "VCFImputable.h"
#include "ParentNoRefImputer.h"
#include "ParentRefImputer.h"
#include "ProgenyRefImputer.h"

class Map;


//////////////////// VCFLargeNoRef ////////////////////

class VCFLargeNoRef : public VCFImputable {
	std::vector<VCFFamilyRecord *>	records;
	ParentNoRefImputer	mat_imputer;
	ParentRefImputer	pat_imputer;
	ProgenyRefImputer	prog_imputer;
	
public:
	VCFLargeNoRef(const STRVEC& samples,
					const std::vector<VCFFamilyRecord *>& rs,
					const Map& map_, double w, const VCFGeno *ref_vcf):
								VCFImputable(samples, ref_vcf->get_ref_vcf()),
								records(rs),
								mat_imputer(records, map_,
											true, ref_vcf, w),
								pat_imputer(records, map_,
											true, ref_vcf, w),
								prog_imputer(records, map_, w) { }
	~VCFLargeNoRef();
	
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
	STRVEC imputed_samples() const override { return samples; }
	void impute() override;
};

#endif
