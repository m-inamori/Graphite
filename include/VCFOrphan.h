#ifndef __VCFORPHAN
#define __VCFORPHAN

#include "VCF.h"

class Map;
class OrphanImputer;


//////////////////// VCFOrphan ////////////////////

class VCFOrphan : public VCFBase, public VCFSmallBase {
private:
	std::vector<VCFRecord *>	records;
	OrphanImputer	*imputer;
	
public:
	VCFOrphan(const std::vector<STRVEC>& header, const STRVEC& s,
						const std::vector<VCFRecord *>& records,
						const std::vector<std::vector<int>>& ref_haps,
						const Map& map_, double w);
	~VCFOrphan();
	
	///// virtual methods for VCFSmallBase /////
	const std::vector<STRVEC>& get_header() const {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const { return VCFBase::get_samples(); }
	std::size_t size() const { return records.size(); }
	VCFRecord *get_record(std::size_t i) const {
		return records[i];
	}
	
	///// non-virtual methods /////
	void impute(std::size_t io);
};
#endif
