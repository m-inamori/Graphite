#include <cmath>
#include "../include/VCFOrphan.h"
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

VCFOrphan::VCFOrphan(const STRVEC& s, const std::vector<GenoRecord *>& rs,
						const std::vector<std::vector<int>>& ref_hs,
						const Map& map_, double w, const VCFSmall *vcf) :
				VCFGenoBase(s, vcf),
				records(rs),
				imputer(records, ref_hs, map_, w) { }

VCFOrphan::~VCFOrphan() {
	Common::delete_all(records);
}

void VCFOrphan::impute_each(size_t i) {
	imputer.impute(i);
}

void VCFOrphan::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	auto	*vcf = c->vcf;
	const size_t	n = vcf->num_samples();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcf->impute_each(i);
	}
}

void VCFOrphan::impute(int T) {
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, this);
	
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
