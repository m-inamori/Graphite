#include <limits>
#include <memory>
#include <cassert>
#include "../include/VCFSelfFillable.h"
#include "../include/VCFSelfHetero.h"
#include "../include/SelfGroups.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSelfFillable ////////////////////

VCFSelfFillable::VCFSelfFillable(const STRVEC& s,
							const std::vector<VCFSelfFillableRecord *>& rs,
							const VCFSmall *vcf) :
										VCFGenoBase(s, vcf), records(rs) { }

VCFSelfFillable::~VCFSelfFillable() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFSelfFillable::modify() {
	const auto	*groups = SelfGroups::create(records);
	const auto	record_sets = groups->create_record_sets();
	for(auto p = record_sets.begin(); p != record_sets.end(); ++p) {
		(*p)->determine_parent_phasing();
	}
	
	delete groups;
	Common::delete_all(record_sets);
}

VCFSelfFillable *VCFSelfFillable::fill(const vector<VCFSelfHetero *>& vcfs,
									const vector<VCFImpSelfRecord *>& records) {
	vector<VCFSelfFillableRecord *>	merged_records
							 = VCFSelfFillable::merge_records(vcfs, records);
	VCFSelfFillable	*vcf = new VCFSelfFillable(vcfs.front()->get_samples(),
												merged_records,
												vcfs.front()->get_ref_vcf());
	vcf->modify();
	return vcf;
}

vector<VCFSelfFillableRecord *> VCFSelfFillable::merge_records(
									const vector<VCFSelfHetero *>& vcfs,
									const vector<VCFImpSelfRecord *>& records) {
	vector<VCFImpSelfRecord *>	all_records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		auto	*vcf = *p;
		auto&	rs = vcf->get_records();
		all_records.insert(all_records.end(), rs.begin(), rs.end());
	}
	all_records.insert(all_records.end(), records.begin(), records.end());
	std::sort(all_records.begin(), all_records.end(),
				[](const VCFImpSelfRecord *r1, const VCFImpSelfRecord *r2) {
					return r1->get_index() < r2->get_index(); });
	const auto&	probs = VCFSelfFillable::calc_probs(all_records, vcfs[0]);
	vector<VCFSelfFillableRecord *>	new_records;
	for(size_t i = 0; i < all_records.size(); ++i) {
		auto	*r = VCFSelfFillableRecord::convert(all_records[i], probs[i]);
		new_records.push_back(r);
	}
	return new_records;
}

vector<vector<VCFRecord::Probs>> VCFSelfFillable::calc_probs(
									const vector<VCFImpSelfRecord *>& records,
									const VCFGenoBase *vcf) {
	const VCFSmall	*orig_vcf = vcf->get_ref_vcf();
	const auto	cols = orig_vcf->extract_columns(vcf->get_samples());
	vector<vector<VCFRecord::Probs>>	prob_table;
	for(size_t i = 0; i < records.size(); ++i) {
		const auto&	geno = records[i]->get_geno();
		const auto	probs = orig_vcf->get_record(i)->parse_PL(geno, cols);
		prob_table.push_back(probs);
	}
	return prob_table;
}
