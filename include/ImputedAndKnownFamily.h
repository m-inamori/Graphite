#ifndef __IMPUTEDANDKNOWNFAMILY
#define __IMPUTEDANDKNOWNFAMILY

#include <array>
#include "VCFImpFamilyRecord.h"
#include "VCFFillableRecord.h"
#include "TypeDeterminer.h"

class VCFOneParentImputed;
class VCFSmallFillable;
class VCFImputable;
class Family;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ImputedAndKnownFamily ////////////////////

namespace ImputedAndKnownFamily {
	std::pair<ParentComb, FillType> classify_record(
											const VCFFamilyRecord *record);
	std::array<std::vector<VCFFillableRecord *>, 4> classify_records(
								const STRVEC& samples,
								const std::vector<VCFFamilyRecord *>& records,
								const VCFSmall *ref_vcf);
	VCFImputable *create(const STRVEC& samples,
									const std::vector<VCFFillableRecord *>& rs,
									bool is_mat_hetero, bool is_mat_imputed,
									const Map& gmap, const VCFSmall *vcf);
	VCFSmallFillable *merge_vcf(const STRVEC& samples,
					const std::array<std::vector<VCFFillableRecord *>, 4>& rss,
					const VCFSmall *vcf);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
							const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	std::size_t compute_upper_NH(const Family *family, std::size_t M,
										std::size_t L, const OptionSmall& op);
	VCFImputable *create_family_vcf(
							const Family *family,
							const std::vector<VCFFamilyRecord *>& records,
							bool is_mat_imputed,
							int num_families,
							const std::vector<std::vector<int>>& ref_haps,
							const VCFSmall *vcf,
							const OptionSmall& op);
	VCFGenoBase *impute_by_parent(
							const VCFSmall *orig_vcf,
							const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const STRVEC& non_imputed_parents,
							const OptionSmall& op);
	inline bool compare_record(const GenoRecord *a, const GenoRecord *b) {
		return a->get_pos() < b->get_pos();
	}
};
#endif
