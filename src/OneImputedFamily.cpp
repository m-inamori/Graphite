#include <algorithm>
#include <cassert>
#include "../include/VCFGeno.h"
#include "../include/OneImputedFamily.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
#include "../include/VCFImputedAndUnknown.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneImputedeFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool OneImputedFamily::is_small(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	N = 1;
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = ((NH*NH << (2*N)) * (2*NH+2*N-1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

bool OneImputedFamily::is_small_ref(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

void OneImputedFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void OneImputedFamily::impute_small_VCFs(vector<VCFImputable *>& v, int T) {
	vector<VCFImputable *>	vcfs(v.begin(), v.end());
	std::sort(vcfs.begin(), vcfs.end(),
				[](const VCFImputable * a, const VCFImputable * b) {
						return a->amount() > b->amount();
				}
	);
	
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

VCFGenoBase *OneImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	L = families.size();
	vector<VCFImputable *>	small_vcfs;
	vector<VCFFamily *>	vcf_garbage;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const bool		is_mat_known = family->is_mat_known();
		const STRVEC&	samples = family->get_samples();
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, 
														orig_vcf, samples);
		if(is_small(ref_haps, L, op)) {
			auto	*vcf1 = new VCFOneParentImputed(samples,
													vcf->get_family_records(),
													ref_haps, is_mat_known,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf1);
		}
		else if(is_small_ref(ref_haps, L, op)) {
			auto	*vcf2 = new VCFOneParentImputedRough(samples,
													vcf->get_family_records(),
													ref_haps, is_mat_known,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf2);
		}
		else {
			// The HMM used here requires little computational amount
			auto	*vcf3 = new VCFImputedAndUnknown(samples,
													vcf->get_family_records(),
													ref_haps, is_mat_known,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf3);
		}
		vcf->clear_records();
		vcf_garbage.push_back(vcf);
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, op.num_threads);
	cout << small_vcfs.size() << " families whose one parent is imputed and"
				<< " the other parent is unknown have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(small_vcfs);
	return new_vcf;
}
