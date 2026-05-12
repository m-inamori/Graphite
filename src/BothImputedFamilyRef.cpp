#include "../include/VCFGeno.h"
#include "../include/BothImputedFamilyRef.h"
#include "../include/VCFBothParentImputed.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/RefCommon.h"
#include "../include/common.h"

using namespace std;


//////////////////// BothImputedFamilyRef ////////////////////

VCFGeno *BothImputedFamilyRef::impute(const VCFSmall *orig_vcf,
									const VCFGenoBase *phased_vcf,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFImputable *>	vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const auto&	samples = family->get_samples();
		VCFGeno	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		const auto	records = RefCommon::merge_family_records(phased_vcf,
																vcf, samples);
		auto	*vcf1 = new VCFBothParentImputed(samples, records, op.map,
													0.01, vcf->get_ref_vcf());
		vcfs.push_back(vcf1);
		delete vcf;
	}
	
	if(vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	cout << vcfs.size()
			<< " both parent imputed families have been imputed." << endl;
	Common::delete_all(vcfs);
	return new_vcf;
}
