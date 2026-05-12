#include <algorithm>
#include <cassert>
#include "../include/OneKnownFamily.h"
#include "../include/VCFOneParentKnown.h"
#include "../include/KnownFamily.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneKnownFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool OneKnownFamily::is_small(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	if(NH >= 100)
		return false;
	
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

// upper NH which passes is_small
size_t OneKnownFamily::compute_upper_NH(const Family *family, size_t M,
										size_t L, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			return NH - 1;
	}
	return 0;	// dummy
}

VCFImputable *OneKnownFamily::create_family_vcf(
									const KnownFamily *family,
									const vector<VCFFamilyRecord *>& records,
									size_t num_families,
									const vector<vector<int>>& ref_haps,
									const VCFSmall *orig_vcf,
									const OptionSmall& op) {
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	const bool	is_mat_known = family->is_mat_known();
	const size_t	M = ref_haps[0].size();
	const size_t	NH = compute_upper_NH(family, M, num_families, op);
	if(is_small(ref_haps, num_families, op)) {
		return new VCFOneParentKnown(family->get_samples(),
										records, ref_haps, is_mat_known,
										op.map, 0.01, orig_vcf);
	}
	else if(NH >= lower_NH) {
		const size_t	NH2 = min(upper_NH, NH);
		const size_t	known_index = is_mat_known ? 0 : 1;
		const auto	gts = VCFFamilyRecord::extract_sample_genotypes(
														known_index, records);
		const auto	ref_haps1 =
				ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH2);
		return new VCFOneParentKnown(family->get_samples(), records,
											ref_haps1, is_mat_known,
											op.map, 0.01, orig_vcf);
	}
	else {
		return NULL;
	}
}

VCFGenoBase *OneKnownFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		auto	*vcf1 = create_family_vcf(family, vcf->get_family_records(),
													N, ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	
	cout << vcfs.size() << " families whose one parent is known and"
				<< " the other parent is unknown have been imputed." << endl;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
