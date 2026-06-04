#include <cmath>
#include <cassert>

#include "../include/VCFFillableRecord.h"
#include "../include/VCFOneParentImputedFast.h"
#include "../include/VCFImpHeteroHomo.h"
#include "../include/VCFHeteroImpHomo.h"
#include "../include/VCFSmallFillable.h"
#include "../include/ClassifyRecord.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFOneParentImputedFast ////////////////////

VCFOneParentImputedFast::~VCFOneParentImputedFast() {
	Common::delete_all(records);
}

STRVEC VCFOneParentImputedFast::imputed_samples() const {
	STRVEC	ss = samples;
	const size_t	index = is_mat_imputed ? 0 : 1;
	ss.erase(ss.begin() + index);
	return ss;
}

size_t VCFOneParentImputedFast::amount() const {
	const size_t	N = num_progenies();
	const size_t	M = records.size();
	return N * M;
}

array<vector<VCFFillableRecord *>, 4>
VCFOneParentImputedFast::classify_records() {
	const auto	cols = get_ref_vcf()->extract_columns(samples);
	// hetero x hetero, homo x hetero, hetero x homo, homo x homo
	array<vector<VCFFillableRecord *>, 4>	rss;
	for(size_t index = 0; index < records.size(); ++index) {
		VCFFamilyRecord	*record = records[index];
		const vector<int>&	geno = record->get_genos();
		const auto	p = ClassifyRecord::classify_family_record(record);
		const ParentComb	comb = p.first;
		const FillType		type = p.second;
		const VCFRecord		*ref_record = get_ref_vcf()->get_record(index);
		const auto			probs = ref_record->parse_PL(geno, cols);
		auto	r = new VCFFillableRecord(record->get_pos(), geno,
												index, type, comb, probs);
		assert(static_cast<int>(type) < 4);
		rss[static_cast<int>(type)].push_back(r);
	}
	return rss;
}

VCFImputable *VCFOneParentImputedFast::create(
										const vector<VCFFillableRecord *>& rs,
										bool is_mat_hetero) const {
	if(is_mat_hetero == is_mat_imputed)
		return new VCFImpHeteroHomo(samples, rs, is_mat_hetero, gmap, vcf);
	else
		return new VCFHeteroImpHomo(samples, rs, is_mat_hetero, gmap, vcf);
}

bool VCFOneParentImputedFast::compare_record(const GenoRecord *a,
											 const GenoRecord *b) {
	return a->get_pos() < b->get_pos();
}

VCFFillable *VCFOneParentImputedFast::merge_vcf(
					const array<vector<VCFFillableRecord *>, 4>& rss) const {
	vector<VCFFillableRecord *>	rs;
	for(int i = 0; i < 4; ++i) {
		for(auto p = rss[i].begin(); p != rss[i].end(); ++p)
			rs.push_back(*p);
	}
	std::sort(rs.begin(), rs.end(), VCFOneParentImputedFast::compare_record);
	return new VCFSmallFillable(samples, rs, vcf);
}


void VCFOneParentImputedFast::impute() {
	const auto	rss = classify_records();
	const auto&	rs_mat = rss[static_cast<int>(FillType::MAT)];
	auto	*mat_vcf = create(rs_mat, true);
	mat_vcf->impute();
	delete mat_vcf;
	const auto&	rs_pat = rss[static_cast<int>(FillType::PAT)];
	auto	*pat_vcf = create(rs_pat, false);
	pat_vcf->impute();
	delete pat_vcf;
	auto	*merged_vcf = merge_vcf(rss);
	merged_vcf->modify(1);
}
