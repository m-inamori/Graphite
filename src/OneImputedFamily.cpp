#include <algorithm>
#include <cassert>
#include "../include/OneImputedFamily.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
#include "../include/KnownFamily.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneImputedeFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool OneImputedFamily::is_small(const vector<vector<int>>& ref_haps, size_t L) {
	const size_t	N = 1;
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = (NH * NH << (2*N)) * (2*NH + 2*N - 1);
	return R*M < 100000000 && R < 100000 && L*R*M < 1000000000;
													// 10^8 & 10^5 & 10^9
}

bool OneImputedFamily::is_small_ref(const vector<vector<int>>& ref_haps,
																	size_t L) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH * NH * (2*NH - 1);
	return R*M < 100000000 && R < 100000 && L*R*M < 1000000000;
													// 10^8 & 10^5 & 10^9
}

void OneImputedFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void OneImputedFamily::impute_small_VCFs(vector<VCFOneParentImputedBase *>& v,
																		int T) {
	vector<VCFOneParentImputedBase *>	vcfs(v.begin(), v.end());
	std::sort(vcfs.begin(), vcfs.end(),
				[](const VCFOneParentImputedBase * a,
				   const VCFOneParentImputedBase * b) {
						return a->amount() > b->amount();
	});
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, vcfs);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&impute_small_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_small_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
}

VCFSmallBase *OneImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const Map& gmap, int num_threads) {
	const size_t	L = families.size();
	vector<VCFOneParentImputedBase *>	small_vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const bool		is_mat_known = family->is_mat_known();
		const STRVEC&	samples = family->get_samples();
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, 
														orig_vcf, samples);
		if(is_small(ref_haps, L)) {
			auto	*vcf1 = new VCFOneParentImputed(vcf->get_header(),
												samples,
												vcf->get_family_records(),
												ref_haps, is_mat_known,
												gmap, 0.01);
			small_vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		else if(is_small_ref(ref_haps, L)) {
			auto	*vcf1 = new VCFOneParentImputedRough(vcf->get_header(),
												samples,
												vcf->get_family_records(),
												ref_haps, is_mat_known,
												gmap, 0.01);
			small_vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		delete vcf;
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, num_threads);
	vector<const VCFSmallBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(small_vcfs);
	return new_vcf;
}
