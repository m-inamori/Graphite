#include <algorithm>
#include <cassert>
#include "../include/VCFGeno.h"
#include "../include/OneImputedFamily.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
#include "../include/VCFImputedAndUnknown.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// OneImputedFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool OneImputedFamily::is_small(const Family *family,
										const vector<vector<int>>& ref_haps,
										size_t L, const OptionSmall& op) {
	const size_t	N = family->num_progenies();
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	if(NH >= 100)
		return false;
	
	const double	R = ((NH*NH << (2*N)) * (2*NH+2*N-1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

bool OneImputedFamily::is_small_ref(const vector<vector<int>>& ref_haps,
											size_t L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	if(NH >= 100)
		return false;
	
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

VCFImputable *OneImputedFamily::create_family_vcf(
									const KnownFamily *family,
									bool is_mat_imputed,
									const vector<VCFFamilyRecord *>& records,
									size_t num_families,
									const vector<vector<int>>& ref_haps,
									const VCFSmall *vcf,
									const OptionSmall& op) {
	const STRVEC&	samples = family->get_samples();
	if(is_small(family, ref_haps, num_families, op))
		return new VCFOneParentImputed(samples, records, ref_haps,
											is_mat_imputed, op.map, 0.01, vcf);
	else if(is_small_ref(ref_haps, num_families, op))
		return new VCFOneParentImputedRough(samples, records, ref_haps,
											is_mat_imputed, op.map, 0.01, vcf);
	else
		return new VCFImputedAndUnknown(samples, records, ref_haps,
											is_mat_imputed, op.map, 0.01, vcf);
}

VCFGenoBase *OneImputedFamily::impute(const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const OptionSmall& op) {
	const size_t	L = families.size();
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const vector<string>&	imputed_samples = imputed_vcf->get_samples();
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const STRVEC&	samples = family->get_samples();
		const bool	is_mat_imputed =
							set_imputed_samples.find(family->get_mat())
											!= set_imputed_samples.end();
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, 
														orig_vcf, samples);
		auto	*vcf1 = create_family_vcf(family, is_mat_imputed,
											vcf->get_family_records(),
											L, ref_haps, orig_vcf, op);
		vcfs.push_back(vcf1);
		vcf->clear_records();
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size() << " families whose one parent is imputed and"
				<< " the other parent is unknown have been imputed." << endl;
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
