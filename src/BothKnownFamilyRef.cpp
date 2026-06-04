#include <algorithm>
#include <cassert>
#include "../include/VCFFamily.h"
#include "../include/BothKnownFamilyRef.h"
#include "../include/BothKnownFamily.h"
#include "../include/VCFBothKnown.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/RefCommon.h"
#include "../include/common.h"

using namespace std;


//////////////////// BothKnownFamilyRef ////////////////////

VCFGeno *BothKnownFamilyRef::impute(const VCFSmall *orig_vcf,
									const VCFGeno *phased_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	N = families.size();
	vector<VCFImputable *>	vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		const auto&	samples = family->get_samples();
		const VCFGeno	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		auto	records = RefCommon::merge_family_records(
													phased_vcf, vcf, samples);
		auto	*vcf1 = BothKnownFamily::create_vcf(family, records, N,
													ref_haps, orig_vcf, op);
		if(vcf1 != NULL) {
			vcfs.push_back(vcf1);
		}
		else {
			Common::delete_all(records);
		}
		delete vcf;
	}
	
	if(vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	cout << vcfs.size()
			<< " both parent known families have been imputed." << endl;
	Common::delete_all(vcfs);
	return new_vcf;
}
