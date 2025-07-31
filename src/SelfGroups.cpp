// SelfGroups.cpp
#include "../include/SelfGroups.h"
#include "../include/VCFSelfFillableRecord.h"

using namespace std;


//////////////////// SelfGroups ////////////////////

VCFSelfFillableRecord *SelfGroups::find_prev_record(size_t i,
												SelfFillType g) const {
	if(i == 0)
		return NULL;
	
	for(size_t j = i-1; j < i; --j) {	// underflow
		const SelfFillType	key = get_group(j).first;
		const auto&	records = get_group(j).second;
		if(key == g)
			return records.back();
	}
	return NULL;
}

VCFSelfFillableRecord *SelfGroups::find_next_record(size_t i,
												SelfFillType g) const {
	for(size_t j = i+1; j < size(); ++j) {
		const SelfFillType	key = get_group(j).first;
		const auto&	records = get_group(j).second;
		if(key == g)
			return records.front();
	}
	return NULL;
}

vector<SelfRecordSet *> SelfGroups::create_record_sets() const {
	vector<SelfRecordSet *>	record_sets;
	for(size_t i = 0; i < size(); ++i) {
		const SelfFillType	key = get_group(i).first;
		const auto&	records = get_group(i).second;
		if(key == SelfFillType::P01 || key == SelfFillType::FILLED)
			continue;
		auto	*prev_record = find_prev_record(i, SelfFillType::P01);
		auto	*next_record = find_next_record(i, SelfFillType::P01);
		for(auto p = records.begin(); p != records.end(); ++p) {
			record_sets.push_back(new SelfRecordSet(*p, prev_record,
														next_record));
		}
	}
	return record_sets;
}

SelfGroups *SelfGroups::create(const vector<VCFSelfFillableRecord *>& records) {
	vector<Group>	groups;
	auto	current_type = records.front()->get_type();
	vector<VCFSelfFillableRecord *>	group(1U, records.front());
	for(auto p = records.begin() + 1; p != records.end(); ++p) {
		auto	*record = *p;
		if(record->get_type() != current_type) {
			groups.push_back(Group(current_type, group));
			group.clear();
			current_type = record->get_type();
		}
		group.push_back(record);
	}
	groups.push_back(Group(current_type, group));
	return new SelfGroups(groups);
}
