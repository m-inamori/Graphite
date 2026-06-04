#include <algorithm>
#include <cassert>
#include "../include/ParentsKnownProgenyImputedFamily.h"
#include "../include/VCFGeno.h"
#include "../include/VCFParentsKnownProgenyImputed.h"
#include "../include/VCFProgenyImputed.h"
#include "../include/KnownFamily.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ParentsKnownProgenyImputedFamily ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool ParentsKnownProgenyImputedFamily::is_small(
											const vector<vector<int>>& ref_haps,
											int L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	if(NH >= 100)
		return false;
	
	// Allow at most one of the four hidden states to change at each transition.
	const double	R = (NH*NH * (8*NH+4)) / op.precision_ratio;
	return R*M < 1e8 && R < 1e5 && L*R*M < 1e9;
}

size_t ParentsKnownProgenyImputedFamily::compute_upper_NH(size_t M, size_t L,
														const OptionSmall& op) {
	size_t	NH = 2;
	while(true) {
		const double	R = (NH*NH * (8*NH+4)) / op.precision_ratio;
		if(!(R*M < 1e8 && R < 1e5 && L*R*M < 1e9))
			break;
		NH += 1;
	}
	return NH - 1;
}

VCFImputable *ParentsKnownProgenyImputedFamily::create_family_vcf(
									const STRVEC& samples,
									const vector<VCFFamilyRecord *>& records,
									int num_families,
									const vector<vector<int>>& ref_haps,
									const VCFSmall *vcf,
									const OptionSmall& op) {
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	const auto	NH = compute_upper_NH(records.size(), num_families, op);
	if(is_small(ref_haps, static_cast<int>(num_families), op)) {
		return new VCFParentsKnownProgenyImputed(samples, records,
													ref_haps, ref_haps,
													op.map, 0.01, vcf);
	}
	else if(NH >= lower_NH) {
		const int	NH2 = min(upper_NH, NH);
		const auto	gts_mat = VCFFamilyRecord::extract_sample_genotypes(
																	0, records);
		const auto	gts_pat = VCFFamilyRecord::extract_sample_genotypes(
																	1, records);
		const auto	ref_haps_mat =
			ReferenceHaplotype::filter_haplotypes(ref_haps, gts_mat, NH2);
		const auto	ref_haps_pat =
			ReferenceHaplotype::filter_haplotypes(ref_haps, gts_pat, NH2);
		return new VCFParentsKnownProgenyImputed(samples, records,
													ref_haps_mat, ref_haps_pat,
													op.map, 0.01, vcf);
	}
	else {
		// If the amount of calculation is so large, impute only a parent
		return new VCFProgenyImputed(samples, records, ref_haps,
												true, op.map, 0.01, vcf);
	}
}

VCFGenoBase *ParentsKnownProgenyImputedFamily::impute(
									const VCFSmall *orig_vcf,
									const VCFGenoBase *imputed_vcf,
									const vector<const KnownFamily *>& families,
									const vector<STRVEC>& imputed_progenies,
									const vector<vector<int>>& ref_haps,
									const OptionSmall& op) {
	const size_t	L = families.size();
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		const STRVEC	samples { family->get_mat(),
								  family->get_pat(),
								  imputed_progenies[i][0] };
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf, 
														orig_vcf, samples);
		auto	*vcf1 = create_family_vcf(samples, vcf->get_family_records(),
													L, ref_haps, orig_vcf, op);
		vcfs.push_back(vcf1);
		vcf->clear_records();
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty())
		return NULL;
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size()
			<< " families whose progeny is imputed have been imputed." << endl;
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
