#include <algorithm>
#include <cassert>
#include "../include/SelfNonImputedFamily.h"
#include "../include/VCFSelfNoImputed.h"
#include "../include/VCFSelfNoImputedRough.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Genotype.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// SelfNonImputedFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool SelfNonImputedFamily::is_small(const Family *family,
										const vector<vector<int>>& ref_haps,
										size_t L, const OptionSmall& op) {
	const size_t	N = family->num_progenies();
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH+2*N-1) << (N*2)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

bool SelfNonImputedFamily::is_small_ref(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

// upper NH which passes is_small_ref
size_t SelfNonImputedFamily::compute_upper_NH(const Family *family, size_t M,
										size_t L, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			return NH - 1;
	}
	return 0;	// dummy
}

VCFSelfImputable *SelfNonImputedFamily::create_family_vcf(const Family *family,
											const vector<GenoRecord *>& records,
											size_t num_families,
											const vector<vector<int>>& ref_haps,
											const VCFSmall *vcf,
											const OptionSmall& op) {
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	vector<string>	samples = family->get_samples();
	samples.erase(samples.begin() + 1);		// Remove one of the two parents
	const size_t	NH = compute_upper_NH(family, vcf->size(),
													num_families, op);
	if(is_small(family, ref_haps, num_families, op)) {
		return new VCFSelfNoImputed(samples, records, ref_haps,
													op.map, 0.01, vcf);
	}
	else if(is_small_ref(ref_haps, (int)num_families, op)) {
		return new VCFSelfNoImputedRough(samples, records, ref_haps,
													op.map, 0.01, vcf);
	}
	else if(NH >= lower_NH) {
		const size_t	NH3 = min(upper_NH, NH);
		const auto	gts = GenoRecord::extract_sample_genotypes(0, records);
		const auto	filtered_ref_haps =
				ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH3);
		return new VCFSelfNoImputedRough(samples, records,
												filtered_ref_haps,
												op.map, 0.01, vcf);
	}
	else {
		return NULL;
	}
}

VCFGeno *SelfNonImputedFamily::impute(
									const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFSelfImputable *>	vcfs;
	vector<VCFGeno *>	vcf_garbage;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		vector<string>	samples = family->get_samples();
		samples.erase(samples.begin() + 1);		// Remove one of the two parents
		auto	*vcf = VCFGeno::create_by_two_vcfs(imputed_vcf,
													orig_vcf, samples);
		auto	*vcf1 = create_family_vcf(family, vcf->get_records(), N,
													ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		
		vcf_garbage.push_back(vcf);
		delete vcf;
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	VCFSelfImputable::impute_VCFs(vcfs, op.num_threads);
	
	cout << vcfs.size() << " self families have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
