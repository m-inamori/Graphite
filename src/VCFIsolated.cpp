#include <algorithm>
#include <cassert>
#include "../include/VCFIsolated.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFIsolated ////////////////////

VCFIsolated::VCFIsolated(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFRecord *> rs, size_t nis, const Map& m) :
										VCFSmall(h, s, rs), VCFMeasurable(m),
										num_imputed_samples(nis) { }

bool VCFIsolated::is_block(const VCFRecord *record,
							const vector<VCFRecord *>& rs) const {
	const double	length = cM(record->pos()) - cM(rs[0]->pos());
	return length < 1.0 || (rs.size() < 10 && length < 10.0);
}

vector<VCFIsolated *> VCFIsolated::divide_by_cM() const {
	// 1cMを超えても10個以内で10cM以内なら塊とみなす
	vector<VCFIsolated *>	rss;
	vector<VCFRecord *>	rs(1, this->records[0]);
	for(auto p = records.begin() + 1; p != records.end(); ++p) {
		VCFRecord	*record = *p;
		if(!this->is_block(record, rs)) {
			rss.push_back(new VCFIsolated(header, samples, rs,
										num_imputed_samples, get_map()));
			rs.clear();
		}
		rs.push_back(record);
	}
	rss.push_back(new VCFIsolated(header, samples, rs,
								num_imputed_samples, get_map()));
	return rss;
}

int VCFIsolated::get_single_gt(const VCFRecord *record, Haplotype hap) const {
	const size_t	i = hap.first;
	const size_t	j = hap.second;
	return record->get_gt(i).c_str()[j*2] - '0';
}

int VCFIsolated::score_each(Haplotype hap1, Haplotype hap2,
									size_t i, VCFRecord *record) const {
	if(record->is_NA(i))
		return 0;
	
	const int	h1 = get_single_gt(record, hap1);
	const int	h2 = get_single_gt(record, hap2);
	const int	int_gt = record->get_int_gt(i);
	return int_gt == h1 + h2 ? 1 : 0;
}

int VCFIsolated::score(Haplotype hap1, Haplotype hap2, size_t i) const {
	int	total_score = 0;
	for(auto p = records.begin(); p != records.end(); ++p)
		total_score += this->score_each(hap1, hap2, i, *p);
	return total_score;
}

vector<VCFIsolated::HaplotypePair>
VCFIsolated::collect_optimal_haplotype_pairs(size_t i) const {
	int	max_score = 0;
	vector<HaplotypePair>	max_combs;
	for(size_t i1 = num_imputed_samples; i1 < samples.size(); ++i1) {
		for(size_t i2 = num_imputed_samples; i2 < samples.size(); ++i2) {
			if(i1 == i2)
				continue;
			for(size_t k = 0; k < 4; ++k) {
				const int	j1 = k >> 1;
				const int	j2 = k & 1;
				const Haplotype	hap1 = make_pair(i1, j1);
				const Haplotype	hap2 = make_pair(i2, j2);
				const int	s = this->score(hap1, hap2, i);
				if(s > max_score) {
					max_score = s;
					max_combs.clear();
					max_combs.push_back(make_pair(hap1, hap2));
				}
				else if(s == max_score) {
					max_combs.push_back(make_pair(hap1, hap2));
				}
			}
		}
	}
	return max_combs;
}

void VCFIsolated::set_haplotype(HaplotypePair hap, size_t i) {
	const Haplotype	hap1 = hap.first;
	const Haplotype	hap2 = hap.second;
	const size_t	i1 = hap1.first;
	const size_t	j1 = hap1.second;
	const size_t	i2 = hap2.first;
	const size_t	j2 = hap2.second;
	for(auto p = records.begin(); p != records.end(); ++p) {
		VCFRecord	*record = *p;
		const char	gt1 = record->get_v()[i1+9].c_str()[j1*2];
		const char	gt2 = record->get_v()[i2+9].c_str()[j2*2];
		char	GT[] = { gt1, '|', gt2, '\0' };
		record->set_GT(i, GT);
	}
}

int VCFIsolated::match_score(HaplotypePair prev_hap, HaplotypePair hap) const {
	const Haplotype	prev_hap1 = prev_hap.first;
	const Haplotype	prev_hap2 = prev_hap.second;
	const Haplotype	hap1 = hap.first;
	const Haplotype	hap2 = hap.second;
	return (prev_hap1 == hap1 ? 1 : 0) + (prev_hap2 == hap2 ? 1 : 0);
}

vector<VCFIsolated::HaplotypePair> VCFIsolated::collect_max_score(
											const vector<HaplotypePair>& combs,
											HaplotypePair prev_hap) const {
	int	max_score = 0;
	vector<HaplotypePair>	max_combs;
	for(auto p = combs.begin(); p != combs.end(); ++p) {
		const int	score = match_score(prev_hap, *p);
		if(score == max_score) {
			max_combs.push_back(*p);
		}
		else if(score > max_score) {
			max_score = score;
			max_combs.clear();
			max_combs.push_back(*p);
		}
	}
	return max_combs;
}

VCFIsolated::HaplotypePair VCFIsolated::impute_cM_each_sample(
											HaplotypePair prev_hap, size_t i) {
	// とりあえず、総当たりにしてみる
	const auto	combs = this->collect_optimal_haplotype_pairs(i);
	
	// 前との一致度が高い組み合わせを集める
	const auto	filtered_combs = collect_max_score(combs, prev_hap);
	
	// 乱数っぽく決める
	const int	j = records[0]->pos() % filtered_combs.size();
	const HaplotypePair	hap = filtered_combs[j];
	this->set_haplotype(hap, i);
	return hap;
}

vector<VCFIsolated::HaplotypePair> VCFIsolated::impute_cM(
									const vector<HaplotypePair>& prev_haps) {
	vector<HaplotypePair>	haps;
	for(size_t i = 0; i < prev_haps.size(); ++i) {
		const auto	prev_hap = prev_haps[i];
		const auto	hap = this->impute_cM_each_sample(prev_hap, i);
		haps.push_back(hap);
	}
	return haps;
}

void VCFIsolated::impute() {
	const Haplotype	h(0, 2);
	vector<HaplotypePair>	haps(num_imputed_samples, make_pair(h, h));
	vector<VCFIsolated *>	vcf_cMs = this->divide_by_cM();
	for(auto p = vcf_cMs.begin(); p != vcf_cMs.end(); ++p) {
		VCFIsolated	*vcf_cM = *p;
		haps = vcf_cM->impute_cM(haps);
		vcf_cM->records.clear();
		delete vcf_cM;
	}
}

VCFSmall *VCFIsolated::extract_isolated_samples() const {
	const STRVEC	isolated_samples(samples.begin(),
									 samples.begin() + num_imputed_samples);
	return extract_samples(isolated_samples);
}

vector<VCFIsolated *> VCFIsolated::create(const VCFSmall *orig_vcf,
											const VCFSmall *imputed_vcf,
											const vector<string>& samples,
											const vector<string>& references,
											const Map& gmap, int num_threads) {
	const auto	sample_columns = orig_vcf->extract_columns(samples);
	const auto	ref_columns = imputed_vcf->extract_columns(references);
	const auto	column_table = divide_columns(sample_columns, num_threads);
	vector<VCFIsolated *>	vcfs;
	for(auto p = column_table.begin(); p != column_table.end(); ++p) {
		const auto	cs = *p;
		STRVEC	new_samples;
		for(auto q = cs.begin(); q != cs.end(); ++q)
			new_samples.push_back(orig_vcf->get_samples()[*q-9]);
		new_samples.insert(new_samples.end(),
							references.begin(), references.end());
		const auto	header = orig_vcf->create_header(new_samples);
		vector<VCFRecord *>	records;
		// samplesをVCFに持たせるため、先にVCFを作って、後でRecordを追加する
		VCFIsolated	*vcf = new VCFIsolated(header, new_samples,
												records, cs.size(), gmap);
		for(size_t i = 0; i < orig_vcf->size(); ++i) {
			VCFRecord	*record = orig_vcf->get_record(i);
			VCFRecord	*imputed_record = imputed_vcf->get_record(i);
			const STRVEC&	v0 = record->get_v();
			STRVEC	v(v0.begin(), v0.begin() + 9);
			for(auto q = cs.begin(); q != cs.end(); ++q)
				v.push_back(record->get_v()[*q]);
			for(auto q = ref_columns.begin(); q != ref_columns.end(); ++q)
				v.push_back(imputed_record->get_v()[*q]);
			VCFRecord	*new_record = new VCFRecord(v, vcf->get_samples());
			vcf->add_record(new_record);
		}
		vcfs.push_back(vcf);
	}
	return vcfs;
}

vector<vector<size_t>> VCFIsolated::divide_columns(
										const vector<size_t>& cs, int num) {
	vector<vector<size_t>>	css;
	if(cs.size() <= (size_t)num) {
		for(auto p = cs.begin(); p != cs.end(); ++p)
			css.push_back(vector<size_t>(1, *p));
	}
	else {
		css.resize(num);
		for(size_t i = 0; i < cs.size(); ++i)
			css[i%num].push_back(cs[i]);
	}
	return css;
}

vector<VCFSmall *> VCFIsolated::impute_all(const vector<VCFIsolated *>& vcfs,
															int num_threads) {
	impute_parellel(vcfs, num_threads);
	
	vector<VCFSmall *>	new_vcfs;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		new_vcfs.push_back((*p)->extract_isolated_samples());
		delete *p;
	}
	return new_vcfs;
}

void VCFIsolated_impute_in_thread(void *config) {
	const auto	*c = (VCFIsolated::ConfigThread *)config;
	for(size_t i = c->first; i < c->size(); i += c->num_threads) {
		auto	vcf = c->vcfs[i];
		vcf->impute();
	}
}

void VCFIsolated::impute_parellel(const vector<VCFIsolated *>& vcfs,
														int num_threads) {
	const int	T = min((int)vcfs.size(), num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(vcfs, (size_t)i, T);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
						(void *(*)(void *))&VCFIsolated_impute_in_thread,
						(void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		VCFIsolated_impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
}
