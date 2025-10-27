#include <limits>
#include <memory>
#include <cassert>
#include "../include/VCFFillable.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/group.h"

#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFFillable ////////////////////

VCFFillable::VCFFillable(const STRVEC& s,
							const std::vector<VCFFillableRecord *>& rs,
							const VCFSmall *vcf) :
									VCFFamilyBase(s, vcf), records(rs) { }

VCFFillable::~VCFFillable() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFFillable::phase_in_thread(void *config) {
	auto	*c = static_cast<const ConfigThreadPhase *>(config);
	const auto&	record_sets = c->record_sets;
	const size_t	n = record_sets.size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		if(record_sets[i]->record->is_fillable_type())
			record_sets[i]->impute(true);
	}
}

void VCFFillable::modify(int T) {
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets();
	
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
}

void VCFFillable::phase_hetero_hetero() {
	// Group records with type 'IMPUTABLE', 'MAT', 'PAT' and 'FIXED'
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets();
	for(auto p = record_sets.begin(); p != record_sets.end(); ++p) {
		(*p)->impute(false);
	}
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		if(record->is_mat_type())
			this->impute_NA_mat(i);
		else if(record->is_pat_type())
			this->impute_NA_pat(i);
		else if(record->is_fillable_type())
			this->impute_others(i);
	}
	Common::delete_all(record_sets);
	delete groups;
}

template<typename Iter>
VCFFillableRecord *VCFFillable::find_neighbor_same_type_record(
							size_t i, size_t k, Iter first, Iter last) const {
	const FillType	type = records[i]->get_type();
	for(auto p = first; p != last; ++p) {
		auto	*record = *p;
		if(record->get_type() == type && record->is_NA(k))
			return record;
	}
	return NULL;
}

VCFFillableRecord *VCFFillable::find_prev_same_type_record(
													size_t i, size_t k) const {
	if(i == 0U)
		return NULL;
	
	return find_neighbor_same_type_record(i, k, records.rend() - i + 1,
																records.rend());
}

VCFFillableRecord *VCFFillable::find_next_same_type_record(
													size_t i, size_t k) const {
	if(i == records.size() - 1)
		return NULL;
	
	const auto	first = records.begin() + i + 1;
	const auto	last = records.end();
	return find_neighbor_same_type_record(i, k, first, last);
}

const RecordSet *VCFFillable::create_recordset(
										size_t i, size_t k, bool is_mat) const {
	auto	*record = records[i];
	auto	*prev_record = find_prev_same_type_record(i, k);
	auto	*next_record = find_next_same_type_record(i, k);
	if(is_mat)
		return new RecordSet(record, prev_record, next_record, NULL, NULL);
	else
		return new RecordSet(record, NULL, NULL, prev_record, next_record);
}

void VCFFillable::impute_NA_mat_each(size_t i, size_t k) {
	const RecordSet	*rs = this->create_recordset(i, k, true);
	rs->impute_NA_mat_each(k);
	delete rs;
}

// impute N/A genotypes in families with heterozygous mat and homozygous pat
void VCFFillable::impute_NA_mat(size_t i) {
	auto	*record = this->records[i];
	for(size_t k = 2; k != get_samples().size(); ++k) {
		if(record->is_NA(k))
			impute_NA_mat_each(i, k);
	}
}

void VCFFillable::impute_NA_pat_each(size_t i, size_t k) {
	const RecordSet	*rs = this->create_recordset(i, k, false);
	rs->impute_NA_pat_each(k);
	delete rs;
}

void VCFFillable::impute_NA_pat(size_t i) {
	auto	*record = this->records[i];
	for(size_t k = 2; k != get_samples().size(); ++k) {
		if(record->is_NA(k))
			impute_NA_pat_each(i, k);
	}
}

pair<int, int> VCFFillable::find_prev_mat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->mat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_mat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->mat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_prev_pat_from(int i, int c) const {
	for(int k = i - 1; k >= 0; --k) {
		const int	from1 = this->records[k]->pat_from(c);
		if(from1 != 0)
			return pair<int, int>(k, from1);
	}
	return pair<int, int>(-1, 0);
}

pair<int, int> VCFFillable::find_next_pat_from(int i, int c) const {
	for(int k = i + 1; k < (int)this->size(); ++k) {
		const int	from2 = this->records[k]->pat_from(c);
		if(from2 != 0)
			return pair<int, int>(k, from2);
	}
	return pair<int, int>(-1, 0);
}

int VCFFillable::select_from(const pair<int, int>& f1,
							 const pair<int, int>& f2, int i) const {
	const int	i1 = f1.first;
	const int	from1 = f1.second;
	const int	i2 = f2.first;
	const int	from2 = f2.second;
	if(from1 == 0 && from2 == 0) {
		// Decide which Haplotype comes from when there is no before or after
		const auto	*r0 = this->records[i];
		return r0->get_pos() % 2 + 1;
	}
	if(from1 == from2)
		return from1;
	else if(from2 == 0)
		return from1;
	else if(from1 == 0)
		return from2;
	else {
		// Determine which Haplotype is coming from by physical distance
		const auto	*r0 = this->records[i];
		const auto	*r1 = this->records[i1];
		const auto	*r2 = this->records[i2];
		if(r0->get_pos() * 2 <= r1->get_pos() + r2->get_pos())
			return from1;
		else
			return from2;
	}
}

int VCFFillable::find_mat_from(int i, int c) const {
	return this->select_from(this->find_prev_mat_from(i, c),
							 this->find_next_mat_from(i, c), i);
}

int VCFFillable::find_pat_from(int i, int c) const {
	return this->select_from(this->find_prev_pat_from(i, c),
							 this->find_next_pat_from(i, c), i);
}

void VCFFillable::impute_others(int i) {
	auto	*record = this->records[i];
	const bool	mat_homo = record->is_homo(0);
	const bool	pat_homo = record->is_homo(1);
	for(size_t k = 2; k != record->num_samples(); ++k) {
		if(record->is_phased(k))
			continue;
		const int	mat_from = mat_homo ? 1 : this->find_mat_from(i, k);
		const int	pat_from = pat_homo ? 1 : this->find_pat_from(i, k);
		record->set_geno(k, record->gt_from_parent(mat_from, pat_from));
	}
}

VCFFillable *VCFFillable::fill(const vector<VCFHeteroHomo *>& vcfs,
				const vector<VCFImpFamilyRecord *>& records, int num_threads) {
	vector<VCFFillableRecord *>	merged_records
							 = VCFFillable::merge_records(vcfs, records, true);
	VCFFillable	*vcf = new VCFFillable(vcfs[0]->get_samples(),
										merged_records, vcfs[0]->get_ref_vcf());
	vcf->modify(num_threads);
	return vcf;
}

vector<VCFFillableRecord *> VCFFillable::merge_records(
									const vector<VCFHeteroHomo *>& vcfs,
									const vector<VCFImpFamilyRecord *>& records,
									bool all_out) {
	vector<VCFImpFamilyRecord *>	all_records = records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		auto	*vcf = *p;
		auto&	hh_records = vcf->get_records();
		for(auto q = hh_records.begin(); q != hh_records.end(); ++q)
			all_records.push_back(*q);
	}
	std::sort(all_records.begin(), all_records.end(),
				[](const VCFImpFamilyRecord *r1, const VCFImpFamilyRecord *r2) {
					return r1->get_index() < r2->get_index(); });
	
	return VCFFillableRecord::convert(all_records, vcfs[0]);
}

vector<vector<VCFFillableRecord *>> VCFFillable::collect_records(
										const vector<VCFFillable *>& vcfs) {
	vector<vector<VCFFillableRecord *>>	rss;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		vector<VCFFillableRecord *>	rs;
		for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
			rs.push_back((*p)->get_fillable_record(i));
		rss.push_back(rs);
	}
	return rss;
}

pair<STRVEC, vector<vector<pair<size_t, size_t>>>>
VCFFillable::integrate_samples(const vector<STRVEC>& sss,
									const STRVEC& orig_samples) {
	// { sample: (index of family, index of inner family) }
	map<string, vector<pair<size_t, size_t>>>	dic;
	for(size_t i = 0; i < sss.size(); ++i) {
		const vector<string>&	samples = sss[i];
		for(size_t j = 0; j < samples.size(); ++j)
			dic[samples[j]].push_back(make_pair(i, j));
	}
	
	STRVEC	new_samples;
	vector<vector<pair<size_t, size_t>>>	pos_samples;
	for(auto p = orig_samples.begin(); p != orig_samples.end(); ++p) {
		const string&	sample = *p;
		if(dic.find(sample) != dic.end()) {
			new_samples.push_back(sample);
			pos_samples.push_back(dic[sample]);
		}
	}
	
	return make_pair(new_samples, pos_samples);
}

// Integrate the VCF so that duplicate samples are one
VCFGeno *VCFFillable::integrate(const VCFSmall *ref_vcf,
								const vector<vector<VCFFillableRecord *>>& rss,
								const vector<STRVEC>& sss,
								const STRVEC& orig_samples) {
	const auto	q = integrate_samples(sss, orig_samples);
	const STRVEC&	samples = q.first;
	// What number of each Family is the sample at?
	const auto&		pos_samples = q.second;
	
	vector<GenoRecord *>	records;
	for(auto p = rss.begin(); p != rss.end(); ++p) {
		auto	*r = VCFFillableRecord::integrate(*p, samples, pos_samples);
		records.push_back(r);
	}
	return new VCFGeno(samples, records, ref_vcf);
}

VCFGeno *VCFFillable::merge(const vector<VCFFillable *>& vcfs,
											const STRVEC& orig_samples) {
	const auto	rss = collect_records(vcfs);
	vector<STRVEC>	sample_table;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		sample_table.push_back((*p)->get_samples());
	}
	return integrate(vcfs[0]->get_ref_vcf(), rss, sample_table, orig_samples);
}

void VCFFillable::fill_in_thread(void *config) {
	const auto	*c = static_cast<const ConfigFillThread *>(config);
	const size_t	n = c->size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	vcfs = c->items[i].first;
		auto	records = c->items[i].second;
		auto	result = VCFFillable::fill(vcfs, records, c->num_threads);
		c->filled_vcfs[i] = result;
	}
}
