#include <algorithm>
#include <cassert>
#include "../include/Orphan.h"
#include "../include/VCFOrphan.h"
#include "../include/Pedigree.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// Orphan ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool Orphan::is_small(const vector<vector<int>>& ref_haps,
												const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R * M < 100000000 && R < 100000;		// 10^8 & 10^5
}

void Orphan::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
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

VCFGenoBase *Orphan::impute(const vector<string>& samples,
								const VCFSmall *orig_vcf,
								const vector<vector<int>>& ref_haps,
								const OptionSmall& op) {
	auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
	if(is_small(ref_haps, op)) {
		auto	*vcf1 = new VCFOrphan(samples, vcf->get_records(),
										ref_haps, op.map, 0.01, orig_vcf);
		impute_small_VCF(vcf1, op.num_threads);
		cout << samples.size() << " orphan samples have been imputed." << endl;
		vcf->clear_records();
		delete vcf;
		return vcf1;
	}
	else {
		delete vcf;
		return NULL;
	}
}
