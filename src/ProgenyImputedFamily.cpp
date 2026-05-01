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

vector<STRVEC> ProgenyImputedFamily::collect_samples(
									const vector<const KnownFamily *>& families,
									const vector<STRVEC>& imputed_progenies) {
	// samples must survive until VCFSmall::join, so create them first.
	const size_t	L = families.size();
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
	return sample_table;
}

VCFGeno *ProgenyImputedFamily::impute(const VCFSmall *orig_vcf,
								const VCFGenoBase *imputed_vcf,
								const vector<const KnownFamily *>& families,
								const vector<vector<string>>& imputed_progenies,
								const vector<vector<int>>& ref_haps,
								const OptionSmall& op) {
	const auto	sample_table = collect_samples(families, imputed_progenies);
	
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const size_t	L = families.size();
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, orig_vcf,
															sample_table[i]);
		auto	*vcf1 = new VCFProgenyImputed(sample_table[i],
											  vcf->get_family_records(),
											  ref_haps, family->is_mat_known(),
											  op.map, 0.01, orig_vcf);
		vcfs.push_back(vcf1);
		vcf->clear_records();
		delete vcf;
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size()
			<< " families whose progeny is imputed have been imputed." << endl;
	
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
