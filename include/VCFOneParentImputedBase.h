#ifndef __ONEIMPUTEDFAMILYBASE
#define __ONEIMPUTEDFAMILYBASE

#include "VCF.h"
#include "VCFFamily.h"

class Family;
class KnownFamily;
class Map;


class VCFOneParentImputedBase : public VCFFamilyBase {
protected:
	const std::vector<VCFFamilyRecord *>	records;
	
public:
	VCFOneParentImputedBase(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							const VCFSmall *vcf) :
									VCFFamilyBase(s, vcf),
									records(rs) { }
	virtual ~VCFOneParentImputedBase() { }
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods /////
	virtual std::size_t amount() const = 0;
	virtual void impute() = 0;
};
#endif
