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

VCFSelfFillable::VCFSelfFillable(const std::vector<STRVEC>& h, const STRVEC& s,
							const std::vector<VCFSelfFillableRecord *>& rs) :
								VCFBase(h, s), VCFSmallBase(), records(rs) { }

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

template<typename Iter>
VCFSelfFillableRecord *VCFSelfFillable::find_neighbor_same_type_record(
							size_t i, size_t c, Iter first, Iter last) const {
	const SelfFillType	type = records[i]->get_type();
	const string&	chromosome = records[i]->chrom();
	for(auto p = first; p != last; ++p) {
		auto	*record = *p;
		if(record->chrom() != chromosome)
			return NULL;
		else if(record->get_type() == type && record->get_GT(c-9) != "./.")
			return record;
	}
	return NULL;
}

VCFSelfFillableRecord *VCFSelfFillable::find_prev_same_type_record(
													size_t i, size_t c) const {
	if(i == 0U)
		return NULL;
	
	return find_neighbor_same_type_record(i, c, records.rend() - i + 1,
																records.rend());
}

VCFSelfFillableRecord *VCFSelfFillable::find_next_same_type_record(
													size_t i, size_t c) const {
	if(i == records.size() - 1)
		return NULL;
	
	const auto	first = records.begin() + i + 1;
	const auto	last = records.end();
	return find_neighbor_same_type_record(i, c, first, last);
}

VCFSelfFillable *VCFSelfFillable::fill(const vector<VCFSelfHetero *>& vcfs,
									const vector<VCFImpSelfRecord *>& records) {
	vector<VCFSelfFillableRecord *>	merged_records
							 = VCFSelfFillable::merge_records(vcfs, records);
	VCFSelfFillable	*vcf = new VCFSelfFillable(vcfs.front()->get_header(),
								vcfs.front()->get_samples(), merged_records);
	vcf->modify();
	return vcf;
}

vector<VCFSelfFillableRecord *> VCFSelfFillable::merge_records(
									const vector<VCFSelfHetero *>& vcfs,
									const vector<VCFImpSelfRecord *>& records) {
	vector<VCFSelfFillableRecord *>	all_records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		auto	*vcf = *p;
		auto&	he_records = vcf->get_records();
		for(auto q = he_records.begin(); q != he_records.end(); ++q)
			all_records.push_back(VCFSelfFillableRecord::convert(*q));
	}
	
	for(auto p = records.begin(); p != records.end(); ++p) {
		all_records.push_back(VCFSelfFillableRecord::convert(*p));
	}
	std::sort(all_records.begin(), all_records.end(),
				[](const VCFSelfFillableRecord *r1, const VCFSelfFillableRecord *r2) {
					return r1->get_index() < r2->get_index(); });
	return all_records;
}

vector<vector<VCFSelfFillableRecord *>> VCFSelfFillable::collect_records(
										const vector<VCFSelfFillable *>& vcfs) {
	vector<vector<VCFSelfFillableRecord *>>	rss;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		vector<VCFSelfFillableRecord *>	rs;
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
			rs.push_back((*p)->get_fillable_record(i));
		rss.push_back(rs);
	}
	return rss;
}
