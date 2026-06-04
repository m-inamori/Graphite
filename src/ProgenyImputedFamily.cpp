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

VCFGeno *ProgenyImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFGenoBase *imputed_vcf,
									const vector<const KnownFamily *>& families,
									const vector<STRVEC>& imputed_progenies,
									const vector<vector<int>>& ref_haps,
									const OptionSmall& op) {
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const size_t	L = families.size();
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		const string&	progeny = imputed_progenies[i][0];
		STRVEC	samples = { family->get_mat(), family->get_pat(), progeny };
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, orig_vcf,
																	samples);
		auto	*vcf1 = new VCFProgenyImputed(samples,
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
	
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
