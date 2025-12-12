#include <algorithm>
#include <cassert>
#include "../include/ImputedAndKnownFamily.h"
#include "../include/VCFImpHeteroHomo.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFOneParentImputedRough.h"
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
		const vector<int>&	geno = record->get_geno();
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

VCFHeteroHomoOnePhased *ImputedAndKnownFamily::create(
										const STRVEC& samples,
										const vector<VCFFillableRecord *>& rs,
										bool is_mat_hetero, bool is_mat_imputed,
										const Map& gmap, const VCFSmall *vcf) {
	if(is_mat_hetero == is_mat_imputed)
		return new VCFImpHeteroHomo(samples, rs, is_mat_hetero, gmap, vcf);
	else
		return new VCFHeteroImpHomo(samples, rs, is_mat_hetero, gmap, vcf);
}

bool ImputedAndKnownFamily::compare_record(const GenoRecord *a,
											const GenoRecord *b) {
	return a->get_pos() < b->get_pos();
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

VCFSmallFillable *ImputedAndKnownFamily::impute(const Family& family,
											VCFFamily *vcf,
											const STRVEC& non_imputed_parents,
											const Map& gmap,
											const VCFSmall *ref_vcf) {
	const STRVEC&	samples = family.get_samples();
	const auto&		records = vcf->get_family_records();
	const auto	rss = classify_records(samples, records, ref_vcf);
	const bool	is_mat_imp = std::find(non_imputed_parents.begin(),
										non_imputed_parents.end(),
										family.get_pat())
												!= non_imputed_parents.end();
	const auto&	rs_mat = rss[static_cast<int>(FillType::MAT)];
	auto	*mat_vcf = create(samples, rs_mat, true, is_mat_imp, gmap, ref_vcf);
	mat_vcf->impute();
	delete mat_vcf;
	const auto&	rs_pat = rss[static_cast<int>(FillType::PAT)];
	auto	*pat_vcf = create(samples, rs_pat, false,
											is_mat_imp, gmap, ref_vcf);
	pat_vcf->impute();
	delete pat_vcf;
	auto	*merged_vcf = merge_vcf(samples, rss, ref_vcf);
	merged_vcf->modify(1);
	return merged_vcf;
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

void ImputedAndKnownFamily::impute_small_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThread *>(config);
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void ImputedAndKnownFamily::impute_small_VCFs(
								vector<VCFImputable *>& vcfs, int T) {
	// VCFOneParentImputed is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	vector<VCFImputable *>	vcfs1(vcfs.begin(), vcfs.end());
	std::sort(vcfs1.begin(), vcfs1.end(),
				[](const VCFImputable * a, const VCFImputable * b) {
					return a->amount() > b->amount();
				}
	);
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, vcfs1);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&impute_small_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_small_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
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
	
	vector<const VCFGenoBase *>	vcfs(N);
	vector<VCFImputable *>	small_vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	// Save ref_haps here as we will impute them all together later.
	vector<vector<vector<int>>>	ref_haps_table(N);
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		const bool	is_mat_imputed = std::find(non_imputed_parents.begin(),
											   non_imputed_parents.end(),
											   family->get_pat())
										!= non_imputed_parents.end();
		const size_t	M = vcf->size();
		const size_t	NH = compute_upper_NH(family, M, N, op);
		if(is_small(family, ref_haps, (int)N, op)) {
			auto	*vcf1 = new VCFOneParentImputed(family->get_samples(),
													vcf->get_family_records(),
													ref_haps, is_mat_imputed,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf1);
			vcfs[i] = vcf1;
			// The records are being reused,
			// so the original VCF is emptied before being deleted.
		}
		else if(is_small_ref(ref_haps, (int)N, op)) {
			auto	*vcf2 = new VCFOneParentImputedRough(family->get_samples(),
													vcf->get_family_records(),
													ref_haps, is_mat_imputed,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf2);
			vcfs[i] = vcf2;
		}
		else if(NH >= lower_NH) {
			const size_t	NH3 = min(upper_NH, NH);
			vector<int>	gts(M);
			size_t	col = is_mat_imputed ? 1 : 0;
			for(size_t j = 0; j < vcf->size(); ++j) {
				gts[j] = vcf->get_record(j)->get_geno()[col];
			}
			// If ref_haps doesn't exist, create it. If it does, reuse it.
			ref_haps_table[i] =
					ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH3);
			auto	*vcf3 = new VCFOneParentImputedRough(vcf->get_samples(),
													vcf->get_family_records(),
													ref_haps_table[i],
													is_mat_imputed,
													op.map, 0.01, orig_vcf);
			small_vcfs.push_back(vcf3);
			vcfs[i] = vcf3;
		}
		else {
			auto	*vcf4 = impute(*family, vcf, non_imputed_parents,
															op.map, orig_vcf);
			vcfs[i] = vcf4;
		}
		vcf->clear_records();
		vcf_garbage.push_back(vcf);
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, op.num_threads);
	cout << N << " families whose one parent is imputed and the other parent"
									<< " is known have been imputed." << endl;
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
