#include <algorithm>
#include <cassert>
#include "../include/VCFGeno.h"
#include "../include/SelfFamilyRef.h"
#include "../include/SelfFamily.h"
#include "../include/VCFSelfImputable.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// SelfFamilyRef ////////////////////

VCFGeno *SelfFamilyRef::impute(const VCFSmall *orig_vcf,
								const VCFGeno *imputed_vcf,
								const vector<vector<int>>& ref_haps,
								const vector<const KnownFamily *>& families,
								const vector<string>& imputed_samples,
								const OptionSmall& op) {
	const size_t	N = families.size();
	if(N == 0)
		return NULL;
	
	vector<VCFSelfImputable *>	vcfs;
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const STRVEC&	samples = family->get_samples();
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		const auto	records = RefCommon::merge_records(imputed_vcf,
															vcf, samples);
		delete vcf;
		auto	*vcf1 = SelfFamily::create_family_vcf(orig_vcf, records,
													  ref_haps, family,
													  set_imputed_samples, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
		}
	}
	
	if(vcfs.empty()) {
		return NULL;
	}
	
	VCFSelfImputable::impute_VCFs(vcfs, op.num_threads);
	
	cout << vcfs.size() << " self families have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
