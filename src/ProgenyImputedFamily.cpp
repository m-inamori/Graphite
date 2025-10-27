#include <algorithm>
#include <cassert>
#include "../include/ProgenyImputedFamily.h"
#include "../include/VCFFamily.h"
#include "../include/VCFProgenyImputed.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ProgenyImputedFamily ////////////////////

void ProgenyImputedFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void ProgenyImputedFamily::impute_small_VCFs(vector<VCFProgenyImputed *>& vcfs,
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

VCFGeno *ProgenyImputedFamily::impute(const VCFSmall *orig_vcf,
								const VCFGenoBase *imputed_vcf,
								const vector<const KnownFamily *>& families,
								const vector<vector<string>>& imputed_progenies,
								const vector<vector<int>>& ref_haps,
								const OptionSmall& op) {
	const size_t	L = families.size();
	// samples must survive until VCFSmall::join, so create them first.
	vector<vector<string>>	sample_table(L);
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		const string	parent = family->is_mat_known() ? family->get_mat() :
														  family->get_pat();
		const string	progeny = imputed_progenies[i][0];
		vector<string>	samples(2);
		samples[0] = parent;
		samples[1] = progeny;
		sample_table[i] = samples;
	}
	
	vector<VCFProgenyImputed *>	small_vcfs;
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, orig_vcf,
															sample_table[i]);
		auto	*vcf1 = new VCFProgenyImputed(sample_table[i],
											  vcf->get_family_records(),
											  ref_haps, family->is_mat_known(),
											  op.map, 0.01, orig_vcf);
		small_vcfs.push_back(vcf1);
		vcf->clear_records();
		delete vcf;
	}
	
	if(small_vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, op.num_threads);
	cout << small_vcfs.size()
			<< " families whose progeny is imputed have been imputed." << endl;
	
	vector<const VCFGenoBase *>	vcfs(small_vcfs.begin(), small_vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(small_vcfs);
	return new_vcf;
}
