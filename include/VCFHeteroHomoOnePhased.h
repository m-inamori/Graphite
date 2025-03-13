#ifndef __VCFHETEROHOMOONEPHASED
#define __VCFHETEROHOMOONEPHASED

#include "VCFImpFamily.h"
#include "VCFFillable.h"
#include "group.h"
#include "option.h"
#include "graph.h"
#include "Map.h"

class Family;
class KnownFamily;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFHeteroHomoOnePhased ////////////////////

class VCFHeteroHomoOnePhased : public VCFBase, public VCFSmallBase,
						public VCFFamilyBase, public VCFMeasurable {
protected:
	std::vector<VCFFillableRecord *>	records;
	const bool	is_mat_hetero;
	
public:
	VCFHeteroHomoOnePhased(const std::vector<STRVEC>& h, const STRVEC& s,
									const std::vector<VCFFillableRecord *>& rs,
									bool mat_hetero, const Map& m) :
												VCFBase(h, s),
												VCFSmallBase(),
												VCFFamilyBase(),
												VCFMeasurable(m),
												records(rs),
												is_mat_hetero(mat_hetero) { }
	virtual ~VCFHeteroHomoOnePhased() { }
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
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
