#ifndef __VCFSELFREF
#define __VCFSELFREF

#include "VCFSelfImputable.h"
#include "ProgenySelfRefImputer.h"
#include "Map.h"

class GenoRecord;
class Map;


//////////////////// VCFSelfRef ////////////////////

class VCFSelfRef : public VCFSelfImputable {
protected:
	std::vector<GenoRecord *>	records;
	ProgenySelfRefImputer	imputer;
	
public:
	VCFSelfRef(const STRVEC& s,
					const std::vector<GenoRecord *>& rs,
					const Map& m, double w, const VCFSmall *vcf);
	VCFSelfRef(const VCFSelfRef&) = delete;
	VCFSelfRef& operator=(const VCFSelfRef&) = delete;
	~VCFSelfRef();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFSelfImputable /////
	std::size_t amount() const override { return 1; }
	void impute() override;
};
#endif
