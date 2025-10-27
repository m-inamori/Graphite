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
class OptionSmall;


//////////////////// ImputedAndKnownFamily ////////////////////

namespace ImputedAndKnownFamily {
	struct ConfigThread {
		std::size_t	first;
		std::size_t	num_threads;
		const std::vector<VCFOneParentImputed *>&	vcfs;
		
		ConfigThread(std::size_t i, std::size_t n,
						const std::vector<VCFOneParentImputed *>& vcfs_) :
									first(i), num_threads(n), vcfs(vcfs_) { }
	};
	
	std::pair<ParentComb, FillType> classify_record(
											const VCFFamilyRecord *record);
	std::array<std::vector<VCFFillableRecord *>, 4> classify_records(
								const STRVEC& samples,
								const std::vector<VCFFamilyRecord *>& records,
								const VCFSmall *ref_vcf);
	VCFHeteroHomoOnePhased *create(const STRVEC& samples,
									const std::vector<VCFFillableRecord *>& rs,
									bool is_mat_hetero, bool is_mat_imputed,
									const Map& gmap, const VCFSmall *vcf);
	bool compare_record(const GenoRecord *a, const GenoRecord *b);
	VCFSmallFillable *merge_vcf(const STRVEC& samples,
					const std::array<std::vector<VCFFillableRecord *>, 4>& rss,
					const VCFSmall *vcf);
	VCFSmallFillable *impute(const Family& family, VCFFamily *vcf,
							const STRVEC& non_imputed_parents,
							const Map& gmap, const VCFSmall *orig_vcf);
	// Is the computational cost sufficiently small even when using ref in HMM?
	bool is_small(const Family *family,
							const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	bool is_small_ref(const std::vector<std::vector<int>>& ref_haps,
							int L, const OptionSmall& op);
	void impute_small_in_thread(void *config);
	void impute_small_VCFs(std::vector<VCFOneParentImputed *>& vcfs, int T);
	VCFGenoBase *impute_by_parent(
							const VCFSmall *orig_vcf,
							const VCFGeno *imputed_vcf,
							const std::vector<std::vector<int>>& ref_haps,
							const std::vector<const KnownFamily *>& families,
							const STRVEC& non_imputed_parents,
							const OptionSmall& op);
};
#endif
