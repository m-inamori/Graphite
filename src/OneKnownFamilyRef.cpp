#include "../include/OneKnownFamilyRef.h"
#include "../include/OneKnownFamily.h"
#include "../include/VCF.h"
#include "../include/VCFImputable.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/RefCommon.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneKnownFamilyRef ////////////////////

VCFGeno *OneKnownFamilyRef::impute(const VCFSmall *orig_vcf,
									const VCFGeno *phased_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFImputable *>	vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const STRVEC&	samples = family->get_samples();
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		const auto	records = RefCommon::merge_family_records(phased_vcf,
																vcf, samples);
		delete vcf;
		auto	*vcf1 = OneKnownFamily::create_family_vcf(family, records,
													N, ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
		}
	}
	
	if(vcfs.empty()) {
		return NULL;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	
	cout << vcfs.size() << " families whose one parent is known and"
				<< " the other parent is unknown have been imputed." << endl;
	
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
