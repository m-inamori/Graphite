#include "../include/OneImputedFamilyRef.h"
#include "../include/OneImputedFamily.h"
#include "../include/VCFImputable.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// OneImputedeFamilyRef ////////////////////

VCFGeno *OneImputedFamilyRef::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	L = families.size();
	vector<VCFImputable *>	vcfs;
	vector<VCFGeno *>	vcf_garbage;
	const vector<string>&	imputed_samples = imputed_vcf->get_samples();
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const STRVEC&	samples = family->get_samples();
		const bool	is_mat_imputed =
							set_imputed_samples.find(family->get_mat())
											!= set_imputed_samples.end();
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		auto	records = RefCommon::merge_family_records(imputed_vcf,
															vcf, samples);
		auto	*vcf1 = OneImputedFamily::create_family_vcf(family,
													is_mat_imputed,
													records, L,
													ref_haps, orig_vcf, op);
		vcfs.push_back(vcf1);
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size() << " families whose one parent is imputed and"
				<< " the other parent is unknown have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
