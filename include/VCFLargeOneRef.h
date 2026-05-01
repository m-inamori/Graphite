#ifndef __VCFLARGEONEREF
#define __VCFLARGEONEREF

#include "VCF.h"
#include "VCFImputable.h"
#include "ParentRefImputer.h"
#include "ProgenyRefImputer.h"

class Map;


//////////////////// VCFLargeOneRef ////////////////////

class VCFLargeOneRef : public VCFImputable {
	std::vector<VCFFamilyRecord *>	records;
	ParentRefImputer	parent_imputer;
	ProgenyRefImputer	prog_imputer;
	
public:
	VCFLargeOneRef(const STRVEC& samples,
					const std::vector<VCFFamilyRecord *>& rs,
					const Map& map_, bool is_mat_ref,
					double w, const VCFGeno *ref_vcf):
								VCFImputable(samples, ref_vcf->get_ref_vcf()),
								records(rs),
								parent_imputer(records, map_,
											is_mat_ref, ref_vcf, w),
								prog_imputer(records, map_, w) { }
	
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
