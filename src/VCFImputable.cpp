#include <algorithm>

#include "../include/VCFImputable.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFImputable ////////////////////

void VCFImputable::impute_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void VCFImputable::impute_VCFs(vector<VCFImputable *>& vcfs, int T) {
	// VCFImputable is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	vector<VCFImputable *>	vcfs1(vcfs.begin(), vcfs.end());
	std::sort(vcfs1.begin(), vcfs1.end(),
				[](const VCFImputable * a, const VCFImputable * b) {
					return a->amount() > b->amount();
				}
	);
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, vcfs1);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&impute_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
}
