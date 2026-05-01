#include <algorithm>

#include "../include/ImputedAndKnownFamilyRef.h"
#include "../include/ImputedAndKnownFamily.h"
#include "../include/VCFImputable.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/RefCommon.h"
#include "../include/common.h"

using namespace std;


//////////////////// ImputedAndKnownFamilyRef ////////////////////

VCFGeno *ImputedAndKnownFamilyRef::impute(
									const VCFSmall *orig_vcf,
									const VCFGeno *phased_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const STRVEC& non_imputed_parents,
									const OptionSmall& op) {
	const size_t	N = families.size();
	if(N == 0)
		return NULL;
	
	vector<VCFImputable *>	vcfs(N);
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const auto&	samples = family->get_samples();
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		auto	records = RefCommon::merge_family_records(phased_vcf,
															vcf, samples);
		const bool	is_mat_imputed = std::find(non_imputed_parents.begin(),
											   non_imputed_parents.end(),
											   family->get_pat())
										!= non_imputed_parents.end();
		vcfs[i] = ImputedAndKnownFamily::create_family_vcf(
												family, records, is_mat_imputed,
												N, ref_haps, orig_vcf, op);
		delete vcf;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << N << " families whose one parent is imputed and the other parent"
									<< " is known have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
