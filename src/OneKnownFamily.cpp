#include <algorithm>
#include <cassert>
#include "../include/OneKnownFamily.h"
#include "../include/VCFOneParentKnown.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneKnownFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool OneKnownFamily::is_small(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

void OneKnownFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute_known_parent();
	}
}

void OneKnownFamily::impute_small_VCFs(vector<VCFOneParentKnown *>& vcfs,
																		int T) {
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

VCFGenoBase *OneKnownFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	L = families.size();
	vector<VCFOneParentKnown *>	vcfs;
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		const bool	is_mat_known = family->is_mat_known();
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		
		if(is_small(ref_haps, L, op)) {
			auto	*vcf1 = new VCFOneParentKnown(family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, is_mat_known,
												   op.map, 0.01, orig_vcf);
			vcf->clear_records();
			vcfs.push_back(vcf1);
		}
		delete vcf;
	}
	
	impute_small_VCFs(vcfs, op.num_threads);
	
	vector<const VCFGenoBase *>	vcfs2;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*vcf1 = *p;
		const string	parent = vcf1->get_known_parent();
		auto	*vcf2 = vcf1->extract_by_samples(vector<string>(1, parent));
		vcfs2.push_back(vcf2);
		delete vcf1;
	}
	
	if(vcfs.empty())
		return NULL;
	
	cout << families.size() << " families whose one parent is known and"
				<< " the other parent is unknown have been imputed." << endl;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	auto	*new_vcf = VCFGeno::join(vcfs2, orig_vcf->get_samples());
	Common::delete_all(vcfs2);
	return new_vcf;
}
