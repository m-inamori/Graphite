#include <algorithm>
#include "../include/LargeFamilyRef.h"
#include "../include/LargeFamily.h"
#include "../include/VCF.h"
#include "../include/VCFFillable.h"
#include "../include/VCFLargeBothRef.h"
#include "../include/VCFLargeNoRef.h"
#include "../include/VCFLargeOneRef.h"
#include "../include/KnownFamily.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// LargeFamilyRef ////////////////////

VCFImputable *LargeFamilyRef::create_family_vcf(const VCFFamilyBase *vcf_family,
														const VCFGeno *ref_vcf,
														const Map& gmap) {
	const auto&	samples = vcf_family->get_samples();
	const auto&	ref_samples = ref_vcf->get_samples();
	auto	records = RefCommon::merge_family_records(ref_vcf,
														vcf_family, samples);
	const bool	is_mat_ref = std::find(ref_samples.begin(), ref_samples.end(),
										vcf_family->mat()) != ref_samples.end();
	const bool	is_pat_ref = std::find(ref_samples.begin(), ref_samples.end(),
										vcf_family->pat()) != ref_samples.end();
	if(is_mat_ref && is_pat_ref)
		return new VCFLargeBothRef(samples, records, gmap,
										0.01, ref_vcf->get_ref_vcf());
	else if(!is_mat_ref && !is_pat_ref)
		return new VCFLargeNoRef(samples, records, gmap, 0.01, ref_vcf);
	else
		return new VCFLargeOneRef(samples, records, gmap,
										is_mat_ref, 0.01, ref_vcf);
}

VCFGeno *LargeFamilyRef::impute(const vector<const KnownFamily *>& families,
								const VCFSmall *orig_vcf,
								const VCFGeno *ref_vcf,
								const Map& gmap, const Option& op) {
	auto	vcfs_family = LargeFamily::impute_all_families(orig_vcf, families,
																	gmap, op);
	if(vcfs_family.empty())
		return NULL;
	
	vector<VCFImputable *>	vcfs;
	for(auto p = vcfs_family.begin(); p != vcfs_family.end(); ++p) {
		vcfs.push_back(create_family_vcf(*p, ref_vcf, gmap));
		delete *p;
	}
	
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size() << " large families have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
