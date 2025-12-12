#include <algorithm>
#include <cassert>
#include "../include/BothKnownFamily.h"
#include "../include/VCFBothKnown.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/ReferenceHaplotype.h"
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

// upper NH which passes is_small
size_t BothKnownFamily::compute_upper_NH(const Family *family, size_t M,
											size_t L, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			return NH - 1;
	}
	return 0;	// dummy
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

VCFGenoBase *BothKnownFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFBothKnown *>	small_vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	// Save ref_haps for mat and pat
	vector<vector<vector<int>>>	ref_haps_table(N*2);
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		const size_t	M = vcf->size();
		const size_t	NH = compute_upper_NH(family, M, N, op);
		if(is_small(ref_haps, (int)N, op)) {
			auto	*vcf1 = new VCFBothKnown(family->get_samples(),
											 vcf->get_family_records(),
											 ref_haps, ref_haps,
											 op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf1);
			vcf->clear_records();
		}
		else if(NH >= lower_NH) {
			const size_t	NH2 = min(upper_NH, NH);
			const auto	gts_mat = vcf->extract_sample_genotypes(0);
			const auto	gts_pat = vcf->extract_sample_genotypes(1);
			ref_haps_table[i*2] =
				ReferenceHaplotype::filter_haplotypes(ref_haps, gts_mat, NH2);
			ref_haps_table[i*2+1] =
				ReferenceHaplotype::filter_haplotypes(ref_haps, gts_pat, NH2);
			auto	*vcf2 = new VCFBothKnown(family->get_samples(),
											 vcf->get_family_records(),
											 ref_haps_table[i*2],
											 ref_haps_table[i*2+1],
											 op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf2);
			vcf->clear_records();
		}
		vcf_garbage.push_back(vcf);
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, op.num_threads);
	cout << small_vcfs.size()
			<< " families whose parents are known have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
