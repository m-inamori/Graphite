#include <algorithm>
#include <cassert>
#include "../include/SelfNonImputedFamilyRef.h"
#include "../include/SelfNonImputedFamily.h"
#include "../include/VCFSelfNoImputed.h"
#include "../include/VCFSelfNoImputedRough.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Genotype.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// SelfNonImputedFamilyRef ////////////////////

VCFGeno *SelfNonImputedFamilyRef::impute(
									const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFSelfImputable *>	vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		vector<string>	samples = family->get_samples();
		samples.erase(samples.begin() + 1);		// Remove one of the two parents
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		const auto	records = RefCommon::merge_records(imputed_vcf,
															vcf, samples);
		delete vcf;
		auto	*vcf1 = SelfNonImputedFamily::create_family_vcf(family,
													records, families.size(),
													ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
		}
	}
	
	if(vcfs.empty())
		return NULL;
	
	VCFSelfImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size() << " self families have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
