#include <algorithm>
#include <cassert>
#include "../include/ParentProgenyImputedFamily.h"
#include "../include/VCFFamily.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ParentProgenyImputedFamily ////////////////////

VCFGeno *ParentProgenyImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFGenoBase *imputed_vcf,
									const vector<const KnownFamily *>& families,
									const STRVEC& non_imputed_parents,
									const vector<vector<int>>& ref_haps,
									const OptionSmall& op) {
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	set<string>	set_non_imputed(non_imputed_parents.begin(),
								non_imputed_parents.end());
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const auto&	samples = family->get_samples();
		auto	*vcf_family = VCFFamily::create_by_two_vcfs(imputed_vcf,
															orig_vcf, samples);
		const bool	should_impute_mat =
						set_non_imputed.find(family->get_mat()) !=
						set_non_imputed.end();
		
		auto	*vcf1 = new VCFOneParentImputed(samples,
											  vcf_family->get_family_records(),
											  ref_haps, should_impute_mat,
											  op.map, 0.01, orig_vcf);
		vcfs.push_back(vcf1);
		vcf_family->clear_records();
		delete vcf_family;
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size()
			<< " families whose one parent is imputed and"
			<<" one progeny is imputed have been imputed." << endl;
	
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
