#include <algorithm>
#include <cassert>
#include "../include/BothKnownFamily.h"
#include "../include/VCFBothKnown.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
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

void BothKnownFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void BothKnownFamily::impute_small_VCFs(vector<VCFBothKnown *>& v, int T) {
	// VCFBothKnown is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	// The order of the VCFs will affect the results,
	// so copy the VCFs and then sort them.
	vector<VCFBothKnown *>	vcfs(v.begin(), v.end());
	std::sort(vcfs.begin(), vcfs.end(),
				[](const VCFBothKnown * a, const VCFBothKnown * b) {
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

VCFSmallBase *BothKnownFamily::impute(const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFBothKnown *>	small_vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		if(is_small(ref_haps, (int)N, op)) {
			auto	*vcf1 = new VCFBothKnown(vcf->get_header(),
												   family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, op.map, 0.01);
			small_vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		delete vcf;
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, op.num_threads);
	cout << small_vcfs.size()
			<< " families whose parents are known have been imputed." << endl;
	vector<const VCFSmallBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(small_vcfs);
	return new_vcf;
}
