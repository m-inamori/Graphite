#ifndef __ONEIMPUTEDFAMILYBASE
#define __ONEIMPUTEDFAMILYBASE

#include "VCF.h"
#include "VCFFamily.h"

class Family;
class KnownFamily;
class Map;


class VCFOneParentImputedBase : public VCFBase,
								public VCFSmallBase, public VCFFamilyBase {
protected:
	const std::vector<VCFFamilyRecord *>	records;
	
public:
	VCFOneParentImputedBase(const std::vector<STRVEC>& header, const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& records);
	virtual ~VCFOneParentImputedBase() { }
	
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
	
	///// virtual methods /////
	virtual void impute() = 0;
	virtual size_t amount() const = 0;
};
#endif
