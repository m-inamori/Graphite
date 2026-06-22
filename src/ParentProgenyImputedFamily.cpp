#include <algorithm>
#include <cassert>
#include "../include/ParentProgenyImputedFamily.h"
#include "../include/VCFFamily.h"
#include "../include/VCFOneParentProgenyImputed.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ParentProgenyImputedFamily ////////////////////

bool ParentProgenyImputedFamily::is_small(const vector<vector<int>>& ref_haps,
									size_t L, size_t P, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * (NH+P*2+1) << (P*2+1)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

size_t ParentProgenyImputedFamily::compute_upper_num_progenies(
									const vector<vector<int>>& ref_haps,
									size_t L, size_t P, const OptionSmall& op) {
	for(size_t p = P; p > 0; --p) {
		if(is_small(ref_haps, L, p, op))
			return p;
	}
	
	return 0;	// dummy
}

const KnownFamily *ParentProgenyImputedFamily::reduce_progenies(
											size_t p, const KnownFamily *f) {
	vector<const Progeny *>	progs(p);
	for(size_t i = 0; i < p; ++i) {
		progs[i] = f->get_progenies()[i]->copy();
	}
	return new KnownFamily(f->get_mat(), f->get_pat(),
							f->is_mat_known(), f->is_pat_known(), progs);
}

VCFImputable *ParentProgenyImputedFamily::create_family_vcf(
							const KnownFamily *family,
							const vector<VCFFamilyRecord *>& records,
							size_t num_families,
							const vector<vector<int>>& ref_haps,
							bool should_impute_mat, const VCFSmall *orig_vcf,
							const OptionSmall& op) {
	auto	*vcf = new VCFFamily(family->get_samples(), records, orig_vcf);
	const size_t	num_progs = family->num_progenies();
	const size_t	p = compute_upper_num_progenies(ref_haps, num_families,
																num_progs, op);
	if(p < num_progs) {
		const auto	*f = reduce_progenies(p, family);
		auto	*vcf1 = vcf->extract(f->get_samples());
		auto	*new_vcf = new VCFOneParentProgenyImputed(f->get_samples(),
													vcf1->get_family_records(),
													ref_haps, should_impute_mat,
													op.map, 0.01, orig_vcf);
		delete vcf;
		vcf1->clear_records();
		delete vcf1;
		delete f;
		return new_vcf;
	}
	else {
		auto	*new_vcf = new VCFOneParentProgenyImputed(family->get_samples(),
													vcf->get_family_records(),
													ref_haps, should_impute_mat,
													op.map, 0.01, orig_vcf);
		vcf->clear_records();
		delete vcf;
		return new_vcf;
	}
}

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
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
														orig_vcf, samples);
		const bool	should_impute_mat =
						set_non_imputed.find(family->get_mat()) !=
						set_non_imputed.end();
		
		auto	*vcf1 = create_family_vcf(family,
											vcf->get_family_records(),
											families.size(),
											ref_haps, should_impute_mat,
											orig_vcf, op);
		vcfs.push_back(vcf1);
		vcf->clear_records();
		vcf_garbage.push_back(vcf);
	}
	
	if(vcfs.empty()) {
		return NULL;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size()
			<< " families whose one parent and"
			<< " one progeny is imputed have been imputed." << endl;
	
	auto	*new_vcf = VCFImputable::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
