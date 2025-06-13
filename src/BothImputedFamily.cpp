#include <algorithm>
#include <cassert>
#include "../include/BothImputedFamily.h"
#include "../include/VCFBothParentImputed.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// BothImputedFamily ////////////////////

void BothImputedFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void BothImputedFamily::impute_VCFs(vector<VCFBothParentImputed *>& v, int T) {
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, v);
	
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

VCFSmallBase *BothImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFBothParentImputed *>	vcfs1;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		auto	*vcf1 = new VCFBothParentImputed(vcf->get_header(),
												   family->get_samples(),
												   vcf->get_family_records(),
												   op.map, 0.01);
		vcfs1.push_back(vcf1);
		vcf->clear_records();
		delete vcf;
	}
	
	if(vcfs1.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_VCFs(vcfs1, op.num_threads);
	vector<const VCFSmallBase *>	vcfs(vcfs1.begin(), vcfs1.end());
	auto	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	cout << vcfs.size()
			<< " both parent imputed families have been imputed." << endl;
	Common::delete_all(vcfs);
	return new_vcf;
}
