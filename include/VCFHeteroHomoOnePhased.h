#ifndef __VCFHETEROHOMOONEPHASED
#define __VCFHETEROHOMOONEPHASED

#include "VCFGeno.h"
#include "VCFImpFamilyRecord.h"
#include "VCFFillable.h"
#include "group.h"
#include "option.h"
#include "graph.h"
#include "Map.h"


//////////////////// VCFHeteroHomoOnePhased ////////////////////

class VCFHeteroHomoOnePhased : public VCFFamilyBase, public VCFMeasurable {
protected:
	std::vector<VCFFillableRecord *>	records;
	const bool	is_mat_hetero;
	
public:
	VCFHeteroHomoOnePhased(const STRVEC& s,
							const std::vector<VCFFillableRecord *>& rs,
							bool mat_hetero, const Map& m,
							const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf),
										VCFMeasurable(m),
										records(rs),
										is_mat_hetero(mat_hetero) { }
	virtual ~VCFHeteroHomoOnePhased() { }
	
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
	const std::vector<VCFFillableRecord *>& get_records() const {
		return records;
	}
	
	virtual void impute() = 0;
};
#endif
