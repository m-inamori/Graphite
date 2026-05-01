#include <algorithm>
#include <cassert>
#include "../include/BothKnownFamily.h"
#include "../include/VCFBothKnown.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// BothKnownFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool BothKnownFamily::is_small(const vector<vector<int>>& ref_haps,
												int L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

// upper NH which passes is_small
size_t BothKnownFamily::compute_upper_NH(const Family *family, size_t M,
											size_t L, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			return NH - 1;
	}
	return 0;	// dummy
}

VCFBothKnown *BothKnownFamily::create_vcf(
								const KnownFamily *family,
								const vector<VCFFamilyRecord *>& records,
								size_t num_families,
								const vector<vector<int>>& ref_haps,
								const VCFSmall *orig_vcf,
								const OptionSmall &op) {
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	const size_t	M = records.size();
	const size_t	NH = compute_upper_NH(family, M, num_families, op);
	if(is_small(ref_haps, (int)num_families, op)) {
		return new VCFBothKnown(family->get_samples(), records,
											 ref_haps, ref_haps,
											 op.map, 0.01, orig_vcf);
	}
	else if(NH >= lower_NH) {
		const size_t	NH2 = min(upper_NH, NH);
		const auto	gts_mat = VCFFamilyRecord::extract_sample_genotypes(
																	0, records);
		const auto	gts_pat = VCFFamilyRecord::extract_sample_genotypes(
																	1, records);
		const auto	ref_haps_mat =
			ReferenceHaplotype::filter_haplotypes(ref_haps, gts_mat, NH2);
		const auto	ref_haps_pat =
			ReferenceHaplotype::filter_haplotypes(ref_haps, gts_pat, NH2);
		return new VCFBothKnown(family->get_samples(), records,
										 ref_haps_mat, ref_haps_pat,
										 op.map, 0.01, orig_vcf);
	}
	else {
		return NULL;
	}
}

VCFGenoBase *BothKnownFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFImputable *>	small_vcfs;
	vector<VCFFamily *>	vcf_garbage;
	// Save ref_haps for mat and pat
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		auto	*vcf1 = create_vcf(family, vcf->get_family_records(),
												N, ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			small_vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		vcf_garbage.push_back(vcf);
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(small_vcfs, op.num_threads);
	cout << small_vcfs.size()
			<< " families whose parents are known have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
