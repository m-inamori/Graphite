#include <map>
#include <set>
#include <cassert>
#include "../include/VCFImpSelfRecord.h"
#include "../include/common.h"

using namespace std;

int VCFImpSelfRecord::get_records_type(const vector<RecordWithPos>& records) {
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

size_t VCFImpSelfRecord::find_fixed_index(
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
