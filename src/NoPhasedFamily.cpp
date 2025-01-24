#include <algorithm>
#include <cassert>
#include "../include/NoPhasedFamily.h"
#include "../include/VCFNoParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
#include "../include/pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/common.h"

using namespace std;


//////////////////// NoPhasedFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool NoPhasedFamily::is_small(const vector<vector<int>>& ref_haps) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH * NH * (2*NH - 1);
	return R * M < 100000000 && R < 100000;		// 10^8 & 10^5
}

VCFSmallBase *NoPhasedFamily::impute(
									const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const Map& gmap) {
	const size_t	N = families.size();
	vector<const VCFSmallBase *>	vcfs;
//	vector<VCFNoParentImputed *>	small_vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		if(is_small(ref_haps)) {
			auto	*vcf1 = new VCFNoParentImputed(vcf->get_header(),
												   family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, gmap, 0.01);
			vcf1->impute();
			auto	*vcf2 = new VCFOneParentImputedRough(vcf->get_header(),
												   family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, true, gmap, 0.01);
			vcf1->clear_records();
			delete vcf1;
			vcf2->impute();
//			small_vcfs.push_back(vcf2);
			vcfs.push_back(vcf2);
			vcf->clear_records();
		}
		delete vcf;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
//	impute_small_VCFs(small_vcfs, num_threads);
	auto	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
