#include <map>
#include <set>
#include <cassert>
#include "../include/VCFImpFamilyRecord.h"
#include "../include/common.h"

using namespace std;

void VCFImpFamilyRecord::set_00x11_parents(size_t i, int gt) {
	const size_t	j = i == 0 ? 1 : 0;
	geno[i] = gt;		// 2 * 3 + 4
	geno[j] = 7 - gt;	// 2 * 3
	const int	prog_gt = (gt == 0) ^ (i == 0) ? 5 : 6;
	for(size_t k = 2; k < geno.size(); ++k) {
		geno[k] = prog_gt;
	}
}

int VCFImpFamilyRecord::get_records_type(const vector<RecordWithPos>& records) {
	int	type = -1;
	for(auto q = records.begin(); q != records.end(); ++q) {
		auto	*record = q->first;
		if(record->is_fixed()) {
			if(type == -1)
				type = 0;
			else if(type == 1)
				return 2;
		}
		else {
			if(type == -1)
				type = 1;
			else if(type == 0)
				return 2;
		}
	}
	return type;
}

size_t VCFImpFamilyRecord::find_fixed_index(
				const vector<std::pair<int, vector<RecordWithPos>>>& items) {
	size_t	all_fixed_index = string::npos;
	for(auto p = items.begin(); p != items.end(); ++p) {
		// 0: all fixed 1: all not fixed 2: otherwise
		const int	records_type = get_records_type(p->second);
		if(records_type == 2)
			return string::npos;
		else if(records_type == 0) {
			if(all_fixed_index == string::npos)
				all_fixed_index = p - items.begin();
			else
				return string::npos;
		}
	}
	return all_fixed_index;
}

std::pair<int, vector<VCFImpFamilyRecord::RecordWithPos>>
VCFImpFamilyRecord::which_is_fixed(const vector<RecordWithPos>& v) {
	vector<int>	gts;
	for(auto p = v.begin(); p != v.end(); ++p) {
		gts.push_back(p->first->unphased(p->second));
	}
	if(Common::is_all_same(gts))
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	// collect reads in Genotype
	// ignore ./.
	map<int, vector<RecordWithPos>>	dic;
	for(size_t k = 0U; k < v.size(); ++k) {
		dic[gts[k]].push_back(v[k]);
	}
	const vector<std::pair<int, vector<RecordWithPos>>>
	items(dic.begin(), dic.end());
	
	// Is the number of the fixed records one?
	vector<bool>	bs(items.size(), true);
	for(size_t i = 0; i < items.size(); ++i) {
		const auto&	w = items[i].second;
		for(auto p = w.begin(); p != w.end(); ++p) {
			if(!p->first->is_fixed())
				bs[i] = false;
		}
	}
	int	counter = 0;
	for(auto p = bs.begin(); p != bs.end(); ++p) {
		if(*p)
			counter += 1;
	}
	if(counter != 1) {
		return make_pair(gts.front(), vector<RecordWithPos>());
	}
	
	// What number is all fixed?
	const size_t	fixed_index = find_fixed_index(items);
	if(fixed_index == string::npos)
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	const int	fixed_GT = items[fixed_index].first;
	if(fixed_GT == 1)	// 0/1
		return make_pair(gts.front(), vector<RecordWithPos>());
	
	vector<RecordWithPos>	wrongs;
	for(size_t i = 0U; i < items.size(); ++i) {
		if(i == (size_t)fixed_index)
			continue;
		const auto&	records = items[i].second;
		for(auto p = records.begin(); p != records.end(); ++p) {
			if(p->first->is_00x11())
				wrongs.push_back(*p);
		}
	}
	return make_pair(fixed_GT, wrongs);
}

void VCFImpFamilyRecord::modify_00x11(const vector<RecordWithPos>& v) {
	const auto	pair = VCFImpFamilyRecord::which_is_fixed(v);
	const int	right_gt = pair.first;
	const vector<VCFImpFamilyRecord::RecordWithPos>& wrongs = pair.second;
	for(auto p = wrongs.begin(); p != wrongs.end(); ++p) {
		auto	*record = p->first;
		const size_t	i = p->second;
		record->set_00x11_parents(i, right_gt);
	}
}
