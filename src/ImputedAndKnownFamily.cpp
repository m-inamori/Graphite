#include <algorithm>
#include <cassert>
#include "../include/ImputedAndKnownFamily.h"
#include "../include/VCFImpHeteroHomo.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
#include "../include/VCFOneParentImputedFast.h"
#include "../include/VCFHeteroImpHomo.h"
#include "../include/VCFSmallFillable.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Genotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ImputedAndKnownFamily ////////////////////

std::pair<ParentComb, FillType> ImputedAndKnownFamily::classify_record(
												const VCFFamilyRecord *record) {
	if(record->is_NA(0) || record->is_NA(1))
		return make_pair(ParentComb::PNA, FillType::IMPUTABLE);
	
	if(record->is_mat_hetero()) {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P01x01, FillType::IMPUTABLE);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x01, FillType::MAT);
		else
			return make_pair(ParentComb::P01x11, FillType::MAT);
	}
	else if(record->is_00(0)) {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P00x01, FillType::PAT);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x00, FillType::FILLED);
		else
			return make_pair(ParentComb::P00x11, FillType::FILLED);
	}
	else {
		if(record->is_pat_hetero())
			return make_pair(ParentComb::P01x11, FillType::PAT);
		else if(record->is_00(1))
			return make_pair(ParentComb::P00x11, FillType::FILLED);
		else
			return make_pair(ParentComb::P11x11, FillType::FILLED);
	}
}

array<vector<VCFFillableRecord *>, 4>
ImputedAndKnownFamily::classify_records(const STRVEC& samples,
								const vector<VCFFamilyRecord *>& records,
								const VCFSmall *ref_vcf) {
	const auto	cols = ref_vcf->extract_columns(samples);
	// hetero x hetero, homo x hetero, hetero x homo, homo x homo
	array<vector<VCFFillableRecord *>, 4>	rss;
	for(size_t index = 0; index < records.size(); ++index) {
		VCFFamilyRecord	*record = records[index];
		const vector<int>&	geno = record->get_genos();
		const auto	p = classify_record(record);
		const ParentComb	comb = p.first;
		const FillType		type = p.second;
		const VCFRecord		*ref_record = ref_vcf->get_record(index);
		const auto			probs = ref_record->parse_PL(geno, cols);
		auto	r = new VCFFillableRecord(record->get_pos(), geno,
												index, type, comb, probs);
		assert(static_cast<int>(type) < 4);
		rss[static_cast<int>(type)].push_back(r);
	}
	return rss;
}

VCFImputable *ImputedAndKnownFamily::create(
										const STRVEC& samples,
										const vector<VCFFillableRecord *>& rs,
										bool is_mat_hetero, bool is_mat_imputed,
										const Map& gmap, const VCFSmall *vcf) {
	if(is_mat_hetero == is_mat_imputed)
		return new VCFImpHeteroHomo(samples, rs, is_mat_hetero, gmap, vcf);
	else
		return new VCFHeteroImpHomo(samples, rs, is_mat_hetero, gmap, vcf);
}

VCFSmallFillable *ImputedAndKnownFamily::merge_vcf(
							const STRVEC& samples,
							const array<vector<VCFFillableRecord *>, 4>& rss,
							const VCFSmall *vcf) {
	vector<VCFFillableRecord *>	rs;
	for(int i = 0; i < 4; ++i) {
		for(auto p = rss[i].begin(); p != rss[i].end(); ++p)
			rs.push_back(*p);
	}
	std::sort(rs.begin(), rs.end(), ImputedAndKnownFamily::compare_record);
	return new VCFSmallFillable(samples, rs, vcf);
}

// Is the computational cost sufficiently small even when using ref in HMM?
bool ImputedAndKnownFamily::is_small(const Family *family,
										const vector<vector<int>>& ref_haps,
										int L, const OptionSmall& op) {
	const size_t	N = family->num_progenies();
	if(N > 2)
		return false;
	
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH+2*N-1) << (N*2)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

bool ImputedAndKnownFamily::is_small_ref(const vector<vector<int>>& ref_haps,
												int L, const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R * M < 1e8 && R < 1e5 && L * R * M < 1e9;
}

size_t ImputedAndKnownFamily::compute_upper_NH(const Family *family, size_t M,
											size_t L, const OptionSmall& op) {
	size_t	NH;
	for(NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 1e8 && R < 1e5 && L * R * M < 1e9))
			break;
	}
	return NH - 1;
}

VCFImputable *ImputedAndKnownFamily::create_family_vcf(
									const Family *family,
									const vector<VCFFamilyRecord *>& records,
									bool is_mat_imputed,
									int num_families,
									const vector<vector<int>>& ref_haps,
									const VCFSmall *orig_vcf,
									const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	const size_t	NH = compute_upper_NH(family, M, num_families, op);
	if(is_small(family, ref_haps, num_families, op)) {
		return new VCFOneParentImputed(family->get_samples(), records,
												ref_haps, is_mat_imputed,
												op.map, 0.01, orig_vcf);
	}
	else if(is_small_ref(ref_haps, num_families, op)) {
		return new VCFOneParentImputedRough(family->get_samples(), records,
													ref_haps, is_mat_imputed,
													op.map, 0.01, orig_vcf);
	}
	else if(NH >= lower_NH) {
		const size_t	NH3 = min(upper_NH, NH);
		vector<int>	gts(M);
		size_t	col = is_mat_imputed ? 1 : 0;
		for(size_t j = 0; j < M; ++j) {
			gts[j] = records[j]->get_geno(col);
		}
		// If ref_haps doesn't exist, create it. If it does, reuse it.
		const auto	filtered_ref_haps =
					ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH3);
		return new VCFOneParentImputedRough(family->get_samples(),
												records,
												filtered_ref_haps,
												is_mat_imputed,
												op.map, 0.01, orig_vcf);
	}
	else {
		return new VCFOneParentImputedFast(family->get_samples(), records,
											is_mat_imputed, op.map, orig_vcf);
	}
}

VCFGenoBase *ImputedAndKnownFamily::impute_by_parent(
									const VCFSmall *orig_vcf,
									const VCFGeno *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const STRVEC& non_imputed_parents,
									const OptionSmall& op) {
	const size_t	N = families.size();
	if(N == 0)
		return NULL;
	
	vector<VCFImputable *>	vcfs(N);
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		const bool	is_mat_imputed = std::find(non_imputed_parents.begin(),
											   non_imputed_parents.end(),
											   family->get_pat())
										!= non_imputed_parents.end();
		// vcfs[i] reuses the record objects from vcf.
		// Clear the records before deleting vcf;
		// otherwise deleting vcf would also delete the records used by vcfs[i].
		vcfs[i] = create_family_vcf(family, vcf->get_family_records(),
									is_mat_imputed, N, ref_haps, orig_vcf, op);
		vcf->clear_records();
		delete vcf;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << N << " families whose one parent is imputed and the other parent"
									<< " is known have been imputed." << endl;
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
