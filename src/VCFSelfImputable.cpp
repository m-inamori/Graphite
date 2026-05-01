#include <algorithm>

#include "../include/VCFSelfImputable.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfImputable ////////////////////

void VCFSelfImputable::impute_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void VCFSelfImputable::impute_VCFs(vector<VCFSelfImputable *>& vcfs, int T) {
	// VCFSelfImputable is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	vector<VCFSelfImputable *>	vcfs1(vcfs.begin(), vcfs.end());
	std::sort(vcfs1.begin(), vcfs1.end(),
				[](const VCFSelfImputable * a, const VCFSelfImputable * b) {
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
