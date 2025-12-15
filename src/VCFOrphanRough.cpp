#include <cmath>
#include "../include/VCFOrphanRough.h"
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

VCFOrphanRough::VCFOrphanRough(
				const STRVEC& s, const std::vector<GenoRecord *>& rs,
				const std::vector<std::vector<std::vector<int>>>& ref_hs_table,
				const Map& map_, double w, const VCFSmall *vcf) :
							VCFGenoBase(s, vcf),
							records(rs),
							imputers(create_imputers(ref_hs_table, map_, w)) { }

VCFOrphanRough::~VCFOrphanRough() {
	Common::delete_all(records);
	Common::delete_all(imputers);
}

vector<OrphanImputer *> VCFOrphanRough::create_imputers(
							const vector<vector<vector<int>>>& ref_haps_table,
							const Map& map_, double w) {
	vector<OrphanImputer *>	imputers;
	for(auto p = ref_haps_table.begin(); p != ref_haps_table.end(); ++p) {
		imputers.push_back(new OrphanImputer(records, *p, map_, w));
	}
	return imputers;
}

void VCFOrphanRough::impute_each(size_t i) {
	imputers[i]->impute(i);
}

void VCFOrphanRough::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	auto	*vcf = c->vcf;
	const size_t	n = vcf->num_samples();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcf->impute_each(i);
	}
}

void VCFOrphanRough::impute(int T) {
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
