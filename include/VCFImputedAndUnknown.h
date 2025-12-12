#ifndef __VCFIMPUTEDANDUNKNOWN
#define __VCFIMPUTEDANDUNKNOWN

#include "VCFImputable.h"
#include "ProgenyImputerByOneParent.h"
#include "ParentImputerByProgeny.h"
#include "ProgenyImputer.h"

class Map;


//////////////////// VCFImputedAndUnknown ////////////////////

class VCFImputedAndUnknown : public VCFImputable {
private:
	std::vector<VCFFamilyRecord *>	records;
	const std::vector<std::vector<int>>&	ref_haps;
	ProgenyImputerByOneParent	prog_imputer;
	ParentImputerByProgeny		parent_imputer;
	ProgenyImputer				progs_imputer;
	
public:
	VCFImputedAndUnknown(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& records,
							const std::vector<std::vector<int>>& ref_haps,
							bool is_mat_imputed,
							const Map& map_, double w, const VCFSmall *vcf);
	VCFImputedAndUnknown(const VCFImputedAndUnknown&) = delete;
	VCFImputedAndUnknown& operator=(const VCFImputedAndUnknown&) = delete;
	~VCFImputedAndUnknown();
	
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
	std::size_t amount() const override;
	void impute() override;
	
	///// non-virtual methods /////
	void clear_records() { records.clear(); }
};
#endif
