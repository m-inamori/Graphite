#include <algorithm>
#include <cassert>
#include "../include/Orphan.h"
#include "../include/VCFOrphan.h"
#include "../include/Pedigree.h"
#include "../include/common.h"

using namespace std;


//////////////////// Orphan ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool Orphan::is_small(const vector<vector<int>>& ref_haps) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH * NH * (2*NH - 1);
	return R * M < 100000000 && R < 100000;		// 10^8 & 10^5
}

void Orphan::impute_small_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	auto	*vcf = c->vcf;
	const size_t	n = vcf->num_samples();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcf->impute(i);
	}
}

void Orphan::impute_small_VCF(VCFOrphan *vcf, int T) {
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, vcf);
	
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

VCFSmallBase *Orphan::impute(const vector<string>& samples,
								const VCFSmall *orig_vcf,
								const vector<vector<int>>& ref_haps,
								const Map& gmap, int num_threads) {
	auto	*vcf = orig_vcf->select_samples(samples);
	if(is_small(ref_haps)) {
		auto	*vcf1 = new VCFOrphan(vcf->get_header(), samples,
										vcf->get_records(),
										ref_haps, gmap, 0.01);
		impute_small_VCF(vcf1, num_threads);
		vcf->clear_records();
		delete vcf;
		return vcf1;
	}
	else {
		return NULL;
	}
}
