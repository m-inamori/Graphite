#ifndef __ONEPHASEDFAMILY
#define __ONEPHASEDFAMILY

#include <array>
#include "VCFImpFamilyRecord.h"
#include "VCFFillableRecord.h"
#include "TypeDeterminer.h"

class VCFOneParentImputed;
class VCFHeteroHomoOnePhased;
class VCFSmallFillable;
class Family;
class KnownFamily;
class Map;


//////////////////// OnePhasedFamily ////////////////////

namespace OnePhasedFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFOneParentImputed *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFOneParentImputed *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
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
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
					const std::vector<std::vector<int>>& ref_haps);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFOneParentImputed *>& vcfs, int T);
	VCFSmallBase *impute_by_parent(
							const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const STRVEC& non_imputed_parents,
							const Map& gmap, int num_threads);
};
#endif
