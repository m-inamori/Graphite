#include <algorithm>
#include <cassert>
#include "../include/NoPhasedFamily.h"
#include "../include/VCFNoParentImputed.h"
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

void NoPhasedFamily::impute_small_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void NoPhasedFamily::impute_small_VCFs(vector<VCFNoParentImputed *>& v, int T) {
	// VCFNoParentImputed is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	// The order of the VCFs will affect the results,
	// so copy the VCFs and then sort them.
	vector<VCFNoParentImputed *>	vcfs(v.begin(), v.end());
	std::sort(vcfs.begin(), vcfs.end(),
				[](const VCFNoParentImputed * a, const VCFNoParentImputed * b) {
						return a->num_samples() > b->num_samples();
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

VCFSmallBase *NoPhasedFamily::impute(const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const Map& gmap, int num_threads) {
	const size_t	N = families.size();
	vector<VCFNoParentImputed *>	small_vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		if(is_small(ref_haps)) {
			auto	*vcf1 = new VCFNoParentImputed(vcf->get_header(),
												   family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, gmap, 0.01);
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
