#ifndef __ONEPHASEDFAMILY
#define __ONEPHASEDFAMILY

#include <array>
#include "VCFImpFamily.h"
#include "VCFFillableRecord.h"
#include "TypeDeterminer.h"

class VCFHeteroHomoOnePhased;
class VCFSmallFillable;
class Family;
class KnownFamily;
class Map;


//////////////////// OnePhasedFamily ////////////////////

namespace OnePhasedFamily {
	int get_gt_type(const VCFRecord *record, int i);
	std::pair<ParentComb, FillType> classify_record(const VCFRecord *record);
	std::array<std::vector<VCFFillableRecord *>, 4> classify_records(
								const std::vector<VCFFamilyRecord *>& records);
	VCFHeteroHomoOnePhased *create(const std::vector<STRVEC>& header,
									const STRVEC& samples,
									const std::vector<VCFFillableRecord *>& rs,
									bool is_mat_hetero, bool is_mat_imputed,
									const Map& gmap);
	bool compare_record(const VCFRecord *a, const VCFRecord *b);
	VCFSmallFillable *merge_vcf(
					const std::vector<STRVEC>& header, const STRVEC& samples,
					const std::array<std::vector<VCFFillableRecord *>, 4>& rss);
	VCFSmallBase *impute(const Family& family, VCFFamily *vcf,
							const STRVEC& non_imputed_parents, const Map& gmap);
	VCFSmallBase *impute_by_parent(
							const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const std::vector<const KnownFamily *>& families,
							const STRVEC& non_imputed_parents,
							const Map& gmap);
};
#endif
