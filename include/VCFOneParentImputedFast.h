#ifndef __VCFONEPARENTIMPUTEDFAST
#define __VCFONEPARENTIMPUTEDFAST

#include "VCFImputable.h"
#include "ParentProgenyImputer.h"
#include "Map.h"

class VCFFillableRecord;
class VCFFillable;


//////////////////// VCFOneParentImputedFast ////////////////////

class VCFOneParentImputedFast : public VCFImputable {
	const std::vector<VCFFamilyRecord *>	records;
	const bool	is_mat_imputed;
	const Map& gmap;
	
public:
	VCFOneParentImputedFast(const STRVEC& s,
							const std::vector<VCFFamilyRecord *>& rs,
							bool is_mat, const Map& map_, const VCFSmall *vcf) :
													VCFImputable(s, vcf),
													records(rs),
													is_mat_imputed(is_mat),
													gmap(map_) { }
	
	~VCFOneParentImputedFast();
	
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
	
private:
	std::array<std::vector<VCFFillableRecord *>, 4> classify_records();
	VCFImputable *create(const std::vector<VCFFillableRecord *>& rs,
													bool is_mat_hetero) const;
	VCFFillable *merge_vcf(
			const std::array<std::vector<VCFFillableRecord *>, 4>& rss) const;
	
private:
	static bool compare_record(const GenoRecord *a, const GenoRecord *b);
};
#endif
