#include <algorithm>
#include <cassert>
#include "../include/OneKnownFamily.h"
#include "../include/VCFOneParentKnown.h"
#include "../include/KnownFamily.h"
#include "../include/ReferenceHaplotype.h"
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

// upper NH which passes is_small
size_t OneKnownFamily::compute_upper_NH(const Family *family, size_t M,
										size_t L, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			return NH - 1;
	}
	return 0;	// dummy
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
	const size_t	N = families.size();
	vector<VCFOneParentKnown *>	vcfs;
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	// Save ref_haps for mat and pat
	vector<vector<vector<int>>>	ref_haps_table(N);
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const bool	is_mat_known = family->is_mat_known();
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		const size_t	M = vcf->size();
		const size_t	NH = compute_upper_NH(family, M, N, op);
		if(is_small(ref_haps, N, op)) {
			auto	*vcf1 = new VCFOneParentKnown(family->get_samples(),
												   vcf->get_family_records(),
												   ref_haps, is_mat_known,
												   op.map, 0.01, orig_vcf);
			vcf->clear_records();
			vcfs.push_back(vcf1);
		}
		else if(NH >= lower_NH) {
			const size_t	NH2 = min(upper_NH, NH);
			const size_t	known_index = is_mat_known ? 0 : 1;
			const auto	gts = vcf->extract_sample_genotypes(known_index);
			ref_haps_table[i] =
					ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH2);
			auto	*vcf2 = new VCFOneParentKnown(family->get_samples(),
													vcf->get_family_records(),
													ref_haps_table[i],
													is_mat_known,
													op.map, 0.01, orig_vcf);
			vcfs.push_back(vcf2);
			vcf->clear_records();
		}
		delete vcf;
	}
	
	impute_small_VCFs(vcfs, op.num_threads);
	
	vector<const VCFGenoBase *>	imputed_vcfs;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*vcf1 = *p;
		const string	parent = vcf1->get_known_parent();
		auto	*imputed_vcf = vcf1->extract_by_samples(vector<string>(1, parent));
		imputed_vcfs.push_back(imputed_vcf);
		delete vcf1;
	}
	
	if(vcfs.empty())
		return NULL;
	
	cout << vcfs.size() << " families whose one parent is known and"
				<< " the other parent is unknown have been imputed." << endl;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	auto	*new_vcf = VCFGeno::join(imputed_vcfs, orig_vcf->get_samples());
	Common::delete_all(imputed_vcfs);
	return new_vcf;
}
