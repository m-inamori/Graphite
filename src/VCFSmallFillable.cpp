#include <limits>
#include <memory>
#include <cassert>
#include "../include/VCFSmallFillable.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/group.h"

#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFSmallFillable ////////////////////

void VCFSmallFillable::phase_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThreadPhase *>(config);
	const auto&	record_sets = c->record_sets;
	const size_t	n = record_sets.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		record_sets[i]->impute(true);
	}
}

void VCFSmallFillable::modify(int T) {
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets_small();
	
	vector<ConfigThreadPhase *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThreadPhase(i, T, record_sets, this);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
			(void *(*)(void *))&phase_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		phase_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	delete groups;
	Common::delete_all(record_sets);
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
	
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->fill_PGT();
}

const RecordSet *VCFSmallFillable::create_recordset(
										size_t i, size_t c, bool is_mat) const {
	auto	*record = records[i];
	auto	*prev_record = find_prev_same_type_record(i, c);
	auto	*next_record = find_next_same_type_record(i, c);
	if(is_mat)
		return new RecordSetSmall(record, prev_record, next_record, NULL, NULL);
	else
		return new RecordSetSmall(record, NULL, NULL, prev_record, next_record);
}
