#include <algorithm>
#include <cassert>
#include "../include/OnePhasedFamily.h"
#include "../include/VCFImpHeteroHomo.h"
#include "../include/VCFOneParentImputed.h"
#include "../include/VCFHeteroImpHomo.h"
#include "../include/VCFSmallFillable.h"
#include "../include/pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/common.h"

using namespace std;


//////////////////// OnePhasedFamily ////////////////////

int OnePhasedFamily::get_gt_type(const VCFRecord *record, int i) {
	if(record->is_hetero(i))
		return 0;	// 0/1
	else if(record->get_gt(i).c_str()[0] == '0')
		return -1;	// 0/0
	else
		return 1;	// 1/1
}

std::pair<ParentComb, FillType> OnePhasedFamily::classify_record(
													const VCFRecord *record) {
	if(record->is_NA(0) || record->is_NA(1))
		return make_pair(ParentComb::PNA, FillType::IMPUTABLE);
	
	const int	mat_gt_type = get_gt_type(record, 0);
	const int	pat_gt_type = get_gt_type(record, 1);
	if(mat_gt_type == 0) {
		switch(pat_gt_type) {
			case 0:  return make_pair(ParentComb::P01x01, FillType::IMPUTABLE);
			case -1: return make_pair(ParentComb::P00x01, FillType::MAT);
			default: return make_pair(ParentComb::P01x11, FillType::MAT);
		}
	}
	else if(mat_gt_type == -1) {
		switch(pat_gt_type) {
			case 0:  return make_pair(ParentComb::P00x01, FillType::PAT);
			case -1: return make_pair(ParentComb::P00x00, FillType::FILLED);
			default: return make_pair(ParentComb::P00x11, FillType::FILLED);
		}
	}
	else {
		switch(pat_gt_type) {
			case 0:  return make_pair(ParentComb::P01x11, FillType::PAT);
			case -1: return make_pair(ParentComb::P00x11, FillType::FILLED);
			default: return make_pair(ParentComb::P11x11, FillType::FILLED);
		}
	}
}

array<vector<VCFFillableRecord *>, 4>
OnePhasedFamily::classify_records(const vector<VCFFamilyRecord *>& records) {
	// hetero x hetero, homo x hetero, hetero x homo, homo x homo
	array<vector<VCFFillableRecord *>, 4>	rss;
	for(size_t index = 0; index < records.size(); ++index) {
		VCFFamilyRecord	*record = records[index];
		const auto	p = classify_record(record);
		const ParentComb	comb = p.first;
		const FillType		type = p.second;
		auto	r = new VCFFillableRecord(record->get_v(),
											record->get_samples(),
											index, type, comb);
		assert(static_cast<int>(type) < 4);
		rss[static_cast<int>(type)].push_back(r);
	}
	return rss;
}

VCFHeteroHomoOnePhased *OnePhasedFamily::create(const vector<STRVEC>& header,
										const STRVEC& samples,
										const vector<VCFFillableRecord *>& rs,
										bool is_mat_hetero, bool is_mat_imputed,
										const Map& gmap) {
	if(is_mat_hetero == is_mat_imputed)
		return new VCFImpHeteroHomo(header, samples, rs, is_mat_hetero, gmap);
	else
		return new VCFHeteroImpHomo(header, samples, rs, is_mat_hetero, gmap);
}

bool OnePhasedFamily::compare_record(const VCFRecord *a, const VCFRecord *b) {
	return a->pos() < b->pos();
}

VCFSmallFillable *OnePhasedFamily::merge_vcf(
					const std::vector<STRVEC>& header, const STRVEC& samples,
					const array<vector<VCFFillableRecord *>, 4>& rss) {
	vector<VCFFillableRecord *>	rs;
	for(int i = 0; i < 4; ++i) {
		for(auto p = rss[i].begin(); p != rss[i].end(); ++p)
			rs.push_back(*p);
	}
	std::sort(rs.begin(), rs.end(), OnePhasedFamily::compare_record);
	return new VCFSmallFillable(header, samples, rs);
}

VCFSmallBase *OnePhasedFamily::impute(const Family& family, VCFFamily *vcf,
											const STRVEC& non_imputed_parents,
											const Map& gmap) {
	const auto&	header = vcf->get_header();
	const STRVEC&	samples = family.get_samples();
	const auto	rss = classify_records(vcf->get_family_records());
	const bool	is_mat_imp = std::find(non_imputed_parents.begin(),
										non_imputed_parents.end(),
										family.get_pat())
												!= non_imputed_parents.end();
	const auto&	rs_mat = rss[static_cast<int>(FillType::MAT)];
	auto	*mat_vcf = create(header, samples, rs_mat, true, is_mat_imp, gmap);
	mat_vcf->impute();
	delete mat_vcf;
	const auto&	rs_pat = rss[static_cast<int>(FillType::PAT)];
	auto	*pat_vcf = create(header, samples, rs_pat, false, is_mat_imp, gmap);
	pat_vcf->impute();
	delete pat_vcf;
	auto	*merged_vcf = merge_vcf(header, samples, rss);
	merged_vcf->modify(1);
	return merged_vcf;
}

// Is the computational cost sufficiently small even when using ref in HMM?
bool OnePhasedFamily::is_small(const Family *family,
								const vector<vector<int>>& ref_haps) {
	const size_t	N = family->num_progenies();
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const size_t	R = NH * NH * (2*NH + 2*N - 1) << (N*2);
	return R * M < 100000000 && R < 100000;		// 10^8 & 10^5
}

void OnePhasedFamily::impute_small_in_thread(void *config) {
	auto	*c = (ConfigThread *)config;
	const auto&	vcfs = c->vcfs;
	const size_t	n = vcfs.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		vcfs[i]->impute();
	}
}

void OnePhasedFamily::impute_small_VCFs(
						vector<VCFOneParentImputed *>& vcfs, int T) {
	// VCFOneParentImputed is heavy for imputation,
	// so make it multi-threaded and impute in order of processing load.
	std::sort(vcfs.begin(), vcfs.end(),
				[](const VCFOneParentImputed * a,
					const VCFOneParentImputed * b) {
						return a->num_progenies() > b->num_progenies();
	});
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(i, T, vcfs);
	
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

VCFSmallBase *OnePhasedFamily::impute_by_parent(
									const VCFSmall *orig_vcf,
									const VCFSmall *imputed_vcf,
									const vector<vector<int>>& ref_haps,
									const vector<const KnownFamily *>& families,
									const STRVEC& non_imputed_parents,
									const Map& gmap, int num_threads) {
	const size_t	N = families.size();
	vector<const VCFSmallBase *>	vcfs(N);
	vector<VCFOneParentImputed *>	small_vcfs;
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf1 = VCFFamily::create_by_two_vcfs(imputed_vcf,
											orig_vcf, family->get_samples());
		if(is_small(family, ref_haps)) {
			const bool	is_mat_imputed = std::find(non_imputed_parents.begin(),
												   non_imputed_parents.end(),
												   family->get_pat())
											!= non_imputed_parents.end();
			auto	*vcf = new VCFOneParentImputed(vcf1->get_header(),
												   family->get_samples(),
												   vcf1->get_family_records(),
												   ref_haps,
												   is_mat_imputed, gmap, 0.01);
			small_vcfs.push_back(vcf);
			vcfs[i] = vcf;
			// The records are being reused,
			// so the original VCF is emptied before being deleted.
			vcf1->clear_records();
		}
		else {
			auto	*imputed_vcf = impute(*family, vcf1,
											non_imputed_parents, gmap);
			vcfs[i] = imputed_vcf;
		}
		delete vcf1;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	impute_small_VCFs(small_vcfs, num_threads);
	auto	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
