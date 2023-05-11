#include <iostream>
#include <cassert>
#include "graphite.h"
#include "SampleManager.h"
#include "VCFOriginal.h"
#include "VCFHeteroHomo.h"
#include "VCFJunkRecord.h"
#include "VCFHomoHomo.h"
#include "VCFHeteroHeteroLite.h"
#include "VCFFillable.h"
#include "VCFHeteroHomoPP.h"
#include "VCFOneParentPhased.h"
#include "VCFProgenyPhased.h"
#include "option.h"
#include "common.h"

using namespace std;


//////////////////// Materials ////////////////////

Materials::~Materials() {
	delete geno_map;
	Common::delete_all(chr_maps);
}

Materials *Materials::create(const Option *option) {
	const auto	*geno_map = Map::read(option->path_map);
	return new Materials(geno_map);
}


//////////////////// process ////////////////////

// HeteroHomoだけ別にする
// このあとHeteroHomoだけ補完するから
// その他はVCFFillableにした後補完する
pair<HeHoRecords, ImpRecords> classify_records(const VCFSmall *vcf,
							const vector<const Family *>& families, double p) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	HeHoRecords	heho_records;
	ImpRecords	other_records;
	for(auto q = families.begin(); q != families.end(); ++q) {
		const Family	*family = *q;
		const Parents	parents = family->parents();
		const STRVEC&	samples = family->get_samples();
		heho_records[parents];
		// vcf_familyはこのループの最後でdeleteするので
		// samplesはfamilyのものを使う
		auto	*vcf_family = VCFFamily::create(vcf, samples);
		const auto	*td = CR->get_TypeDeterminer(family->num_progenies(), p);
		for(size_t i = 0; i < vcf_family->size(); ++i) {
			VCFFamilyRecord	*record = vcf_family->get_record(i);
			const auto	pair1 = CR->classify(record, td);
			const ParentComb	pc = pair1.first;
			const WrongType	wrong_type = pair1.second;
			const STRVEC&	v = record->get_v();
			if(pc == ParentComb::PNA) {		// 候補が無い
				auto	*r = new VCFJunkRecord(v, samples, i, wrong_type);
				other_records[parents].push_back(r);
			}
			else if(TypeDeterminer::is_homohomo(pc)) {
				auto	*r_ = new VCFHomoHomoRecord(v, samples, i,
														wrong_type, pc);
				auto	*r = r_->impute();
				delete r_;
				other_records[parents].push_back(r);
			}
			else if(TypeDeterminer::is_heterohomo(pc)) {
				auto	*r = new VCFHeteroHomoRecord(v, samples, i,
														wrong_type, pc);
				heho_records[parents].push_back(r);
			}
			else {
				auto	*r = new VCFHeteroHeteroLiteRecord(v, samples, i,
															wrong_type, pc);
				other_records[parents].push_back(r);
			}
		}
		delete vcf_family;
	}
	
	return make_pair(heho_records, other_records);
}

vector<vector<VCFImpFamilyRecord *>>
							sort_records(const HeHoRecords& heho_records,
											const ImpRecords& other_records) {
	int	max_index = 0;
	for(auto p = heho_records.begin(); p != heho_records.end(); ++p) {
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			max_index = std::max(max_index, (*q)->get_index());
	}
	for(auto p = other_records.begin(); p != other_records.end(); ++p) {
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			max_index = std::max(max_index, (*q)->get_index());
	}
	
	vector<vector<VCFImpFamilyRecord *>>	records(max_index + 1);
	for(auto p = heho_records.begin(); p != heho_records.end(); ++p) {
		const auto&	rs1 = p->second;
		for(auto q = rs1.begin(); q != rs1.end(); ++q) {
			if((*q)->is_right())
				records[(*q)->get_index()].push_back(*q);
		}
	}
	for(auto p = other_records.begin(); p != other_records.end(); ++p) {
		const auto&	rs1 = p->second;
		for(auto q = rs1.begin(); q != rs1.end(); ++q) {
			if((*q)->is_homohomo())
				records[(*q)->get_index()].push_back(*q);
		}
	}
	return records;
}

// 0/0 x 1/1のrecordで別の家系の親になっているとき、
// 他のホモ×ホモやヘテロ×ホモとGenotypeが違うとき、修正する
void modify_00x11(const HeHoRecords& heho_records, ImpRecords& other_records) {
	auto	records = sort_records(heho_records, other_records);
	
	for(auto p = records.begin(); p != records.end(); ++p) {
		if(p->size() >= 2)
			VCFImpFamilyRecord::modify_00x11(*p);
	}
}

// FamilyごとにVCFHeteroHomoを作って親ごとに格納する
tuple<VCFHeteroHomo *, VCFHeteroHomo *, vector<VCFHeteroHomoRecord *>>
make_VCFHeteroHomo(const vector<VCFHeteroHomoRecord *>& records,
									const Family *family,
									const VCFSmall *vcf, const Map& geno_map) {
	const STRVEC&	samples = family->get_samples();
	const vector<STRVEC>	header = vcf->create_header(samples);
	vector<VCFHeteroHomoRecord *>	heho_mat_records;
	vector<VCFHeteroHomoRecord *>	heho_pat_records;
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		VCFHeteroHomoRecord	*record = *p;
		if(!record->is_imputable())
			unused_records.push_back(record);
		else if(record->is_mat_hetero())
			heho_mat_records.push_back(record);
		else
			heho_pat_records.push_back(record);
	}
	auto	*vcf_mat = new VCFHeteroHomo(header, samples,
											heho_mat_records, geno_map);
	auto	*vcf_pat = new VCFHeteroHomo(header, samples,
											heho_pat_records, geno_map);
	return make_tuple(vcf_mat, vcf_pat, unused_records);
}

void impute_in_thread(void *config) {
	const auto	*c = (ConfigImpThread *)config;
	const size_t	n = c->size();
	// 家系数が少ない場合、スレッドの内部でも並列に処理したい
	const int	num_inner_threads = (c->option->num_threads + n - 1) / n;
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	result = VCFHeteroHomo::impute_vcfs(c->vcfs_heho[i],
												c->option, num_inner_threads);
		c->imputed_vcfs[i] = result;
	}
}

vector<ImpResult> impute_hetero_homo_parellel(
						const map<string, vector<VCFHeteroHomo *>>& vcfs,
						const Option *option) {
	vector<ImpResult>	results(vcfs.size());
	
	// vectorにしたほうがマルチスレッドにしやすい
	vector<vector<VCFHeteroHomo *>>	vcfss;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		vcfss.push_back(p->second);
	}
	
	const int	T = min((int)vcfs.size(), option->num_threads);
	vector<ConfigImpThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigImpThread(vcfss, option, i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&impute_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return results;
}

pair<map<Parents, vector<VCFHeteroHomo *>>, HeHoRecords>
		impute_hetero_homo_core(const HeHoRecords& records,
								const VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option) {
	// あとでFamilyを回復するために必要
	map<Parents, const Family *>	parents_to_family;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const Parents	parents = p->first;
		const Family	*family = sample_man->get_large_family(parents);
		parents_to_family.insert(make_pair(parents, family));
	}
	
	map<string, vector<VCFHeteroHomo *>>	vcfs;
	HeHoRecords	unused_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const Parents&	parents = p->first;
		const Family	*family = parents_to_family[parents];
		const vector<VCFHeteroHomoRecord *>&	heho_records = p->second;
		const auto	t = make_VCFHeteroHomo(heho_records,
											family, vcf, geno_map);
		VCFHeteroHomo	*vcf_mat = get<0>(t);
		VCFHeteroHomo	*vcf_pat = get<1>(t);
		// Hetero x Homoで使われなかったRecordを再利用する
		const vector<VCFHeteroHomoRecord *>&	unused_records1 = get<2>(t);
		vcfs[family->get_mat()].push_back(vcf_mat);
		vcfs[family->get_pat()].push_back(vcf_pat);
		auto&	urs = unused_records[parents];
		urs.insert(urs.end(), unused_records1.begin(), unused_records1.end());
	}
	
	const auto	w = impute_hetero_homo_parellel(vcfs, option);
	
	map<Parents, vector<VCFHeteroHomo *>>	imputed_vcfs;
	for(auto p = w.begin(); p != w.end(); ++p) {
		const vector<VCFHeteroHomo *>&	imputed_vcfs_heho = p->first;
		const vector<VCFHeteroHomoRecord *>&	unused_records2 = p->second;
		for(auto q = imputed_vcfs_heho.begin();
								q != imputed_vcfs_heho.end(); ++q) {
			VCFHeteroHomo	*vcf = *q;
			imputed_vcfs[vcf->parents()].push_back(vcf);
		}
		
		for(auto q = unused_records2.begin(); q != unused_records2.end(); ++q) {
			VCFHeteroHomoRecord	*record = *q;
			record->enable_modification();
			unused_records[record->parents()].push_back(record);
		}
	}
	return make_pair(imputed_vcfs, unused_records);
}

void join_records(ImpRecords& records, HeHoRecords& unused_records) {
	for(auto q = unused_records.begin(); q != unused_records.end(); ++q) {
		const auto&	parents = q->first;
		const auto&	rs = q->second;
		auto&	records_ = records[parents];
		records_.insert(records_.end(), rs.begin(), rs.end());
		sort(records_.begin(), records_.end(),
				[](const VCFFamilyRecord *r1, const VCFFamilyRecord *r2) {
					return r1->pos() < r2->pos(); });
	}
}

pair<map<Parents, vector<VCFHeteroHomo *>>, ImpRecords>
		impute_hetero_homo(const VCFSmall *orig_vcf, SampleManager *sample_man,
									const Map& geno_map, const Option *option) {
	const auto	pair1 = classify_records(orig_vcf,
							sample_man->get_large_families(), option->ratio);
	const auto	heho_records = pair1.first;
	auto	other_records = pair1.second;
	modify_00x11(heho_records, other_records);
	
	// HeteroHomoとそれ以外を分けて、HeteroHomoを補完する
	// 使われなかったHeteroHomoはvillに回す
	const auto	pair2 = impute_hetero_homo_core(heho_records, orig_vcf,
												sample_man, geno_map, option);
	auto	imputed_vcfs = pair2.first;
	auto	unused_records = pair2.second;
	
	join_records(other_records, unused_records);
	
	return make_pair(imputed_vcfs, other_records);
}

void fill_in_thread(void *config) {
	const auto	*c = (ConfigFillThread *)config;
	const size_t	n = c->size();
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	vcfs = c->items[i].first;
		auto	records = c->items[i].second;
		auto	result = VCFFillable::fill(vcfs, records, c->all_out);
		c->filled_vcfs[i] = result;
	}
}

vector<VCFFillable *> fill_parellel(vector<Item>& items, const Option *op) {
	vector<VCFFillable *>	results(items.size());
	
	const int	T = min((int)items.size(), op->num_threads);
	vector<ConfigFillThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigFillThread(items, true, i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&fill_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		fill_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return results;
}

void delete_items(const vector<Item>& items) {
	for(auto q = items.begin(); q != items.end(); ++q) {
		for(auto r = q->first.begin(); r != q->first.end(); ++r)
			delete *r;
		for(auto r = q->second.begin(); r != q->second.end(); ++r)
			delete *r;
	}
}

vector<VCFFillable *> fill(map<Parents, vector<VCFHeteroHomo *>>& imputed_vcfs,
								ImpRecords& other_records, const Option *op) {
	vector<Item>	items;
	for(auto q = imputed_vcfs.begin(); q != imputed_vcfs.end(); ++q) {
		const Parents&	parents = q->first;
		vector<VCFHeteroHomo *>&	vcfs = q->second;
		items.push_back(make_pair(vcfs, other_records[parents]));
	}
	vector<VCFFillable *>	filled_vcfs = fill_parellel(items, op);
	delete_items(items);
	return filled_vcfs;
}

VCFSmall *fill_and_merge_vcf(
					map<Parents, vector<VCFHeteroHomo *>>& imputed_vcfs,
						ImpRecords& other_records,
						const STRVEC& samples, const Option *option) {
	// Familyごとに残りのRecordをphaseする
	vector<VCFFillable *>	filled_vcfs = fill(imputed_vcfs,
												other_records, option);
	// FamilyごとのVCFを全て統合する
	VCFSmall	*vcf_integrated = VCFFillable::merge(filled_vcfs,
														samples, option);
	for(auto p = filled_vcfs.begin(); p != filled_vcfs.end(); ++p)
		delete *p;
	
	return vcf_integrated;
}

VCFRecord *merge_progeny_records(vector<VCFFillable *>& vcfs,
									size_t i, const STRVEC& samples) {
	const STRVEC&	v1 = vcfs.front()->get_records()[i]->get_v();
	STRVEC	v(v1.begin(), v1.begin() + 9);
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const STRVEC&	v2 = (*p)->get_records()[i]->get_v();
		v.insert(v.end(), v2.begin() + 11, v2.end());
	}
	return new VCFRecord(v, samples);
}

VCFSmall *impute_vcf_by_parents(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				const vector<const Family *>& families, const Map& geno_map) {
	vector<VCFFillable *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		auto	*family_vcf = VCFHeteroHomoPP::impute_by_parents(
												orig_vcf, merged_vcf,
												(*p)->get_samples(), geno_map);
		vcfs.push_back(family_vcf);
	}
	
	STRVEC	samples;	// 子どもだけ集める
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFamily	*vcf = *p;
		STRVEC	ss = vcf->get_samples();
		samples.insert(samples.end(), ss.begin() + 2, ss.end());
	}
	
	// samplesの所有権をnew_vcfに持たせる
	auto	new_header = orig_vcf->create_header(samples);
	vector<VCFRecord *>	empty_records;
	auto	*new_vcf = new VCFSmall(new_header, samples, empty_records);
	const auto&	samples_ = new_vcf->get_samples();
	
	vector<VCFRecord *>	merged_records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		new_vcf->add_record(merge_progeny_records(vcfs, i, samples_));
	}
	
	Common::delete_all(vcfs);
	
	return new_vcf;
}

VCFSmall *impute_vcf_by_parent(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				const vector<const Family *>& families,
				const Map& geno_map, const SampleManager *sample_man) {
	vector<VCFFamily *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		const bool	is_mat_phased = sample_man->is_imputed(family->get_mat());
		const auto&	samples = family->get_samples();
		auto	*family_vcf = VCFOneParentPhased::impute_by_parent(
									orig_vcf, merged_vcf,
									samples, is_mat_phased, geno_map);
		vcfs.push_back(family_vcf);
	}
	
	STRVEC	samples;	// 今回phasingした親と子どもを集める
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFamily	*vcf = *p;
		STRVEC	ss = vcf->get_samples();
		samples.insert(samples.end(), ss.begin() + 2, ss.end());
		if(ss[0] == "0" || ss[1] == "0")
			continue;
		else if(sample_man->is_imputed(ss[0]))
			samples.push_back(ss[1]);
		else
			samples.push_back(ss[0]);
	}
	
	// samplesの所有権をnew_vcfに持たせる
	auto	new_header = orig_vcf->create_header(samples);
	vector<VCFRecord *>	empty_records;
	auto	*new_vcf = new VCFSmall(new_header, samples, empty_records);
	const auto&	samples_ = new_vcf->get_samples();
	
	vector<VCFRecord *>	merged_records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		new_vcf->add_record(VCFOneParentPhased::merge_records(
													vcfs, i, samples_));
	}
	
	Common::delete_all(vcfs);
	
	return new_vcf;
}

VCFSmall *impute_vcf_by_progenies(VCFSmall *orig_vcf, VCFSmall *merged_vcf,
									const vector<const Family *>& families,
									SampleManager *sample_man) {
	vector<VCFProgenyPhased *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		vector<size_t>	ppi;
		const vector<string>&	samples = family->get_samples();
		for(size_t i = 0; i < samples.size(); ++i) {
			if(sample_man->is_imputed(samples[i]))
				ppi.push_back(i);
		}
		auto	*vcf = VCFProgenyPhased::impute_by_progeny(orig_vcf, merged_vcf,
													family->get_samples(), ppi);
		vcfs.push_back(vcf);
	}
	
	vector<VCFSmall *>	vcfs2(vcfs.begin(), vcfs.end());
	VCFSmall	*new_vcf = VCFSmall::join(vcfs2, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}

VCFSmall *impute_vcf_chr(VCFSmall *orig_vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option) {
	cerr << "chr : " << orig_vcf->get_records().front()->chrom() << endl;
	
	const auto&	all_samples = orig_vcf->get_samples();
	auto	p = impute_hetero_homo(orig_vcf, sample_man, geno_map, option);
	auto	imputed_vcfs = p.first;
	auto	other_records = p.second;
	auto	*merged_vcf = fill_and_merge_vcf(imputed_vcfs, other_records,
														all_samples, option);
	if(option->only_large_families)
		return merged_vcf;
	
	// 両親が補完されているが子どもが少ない家系を補完する
	// 補完できる家系がなくなるまで繰り返す
	sample_man->add_imputed_samples(merged_vcf->get_samples());
	while(true) {
		// 両親が補完されているが子どもが少ない家系を補完する
		auto	families1 = sample_man->extract_small_families();
		if(!families1.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_parents(orig_vcf,
												merged_vcf, families1, geno_map);
			Common::delete_all(families1);
			vector<VCFSmall *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			continue;
		}
		
		// 片親が補完されている家系を補完する
		auto	families2 = sample_man->extract_single_parent_phased_families();
		if(!families2.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_parent(
											orig_vcf, merged_vcf,
											families2, geno_map, sample_man);
			vector<VCFSmall *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			Common::delete_all(families2);
			continue;
		}
		
		// Impute families whose progenies have been imputed
		auto	families3 = sample_man->extract_progenies_phased_families();
		if(!families3.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_progenies(
														orig_vcf, merged_vcf,
														families3, sample_man);
			vector<VCFSmall *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			Common::delete_all(families3);
			continue;
		}
		else
			break;
	}
	sample_man->clear_imputed_samples();
	return merged_vcf;
}

void impute_VCF(const Option *option) {
	// chromosomeごとに処理する
	Materials	*materials = Materials::create(option);
	VCFHuge	*vcf = VCFHuge::read(option->path_vcf);
	SampleManager	*sample_man = SampleManager::create(
										option->path_ped, vcf->get_samples(),
										option->lower_progs, option->families);
	VCFHuge::ChromDivisor	divisor(vcf);
	bool	first = true;
	for(int	chrom_index = 0; ; ++chrom_index) {
		VCFSmall	*vcf_chrom = divisor.next();
		if(vcf_chrom == NULL)
			break;
		else if(!option->is_efficient_chrom(chrom_index)) {
			delete vcf_chrom;
			continue;
		}
		
		const Map	*gmap = materials->get_chr_map(chrom_index);
		const VCFSmall	*vcf_imputed = impute_vcf_chr(vcf_chrom, sample_man,
																*gmap, option);
		delete vcf_chrom;
		if(first) {
			ofstream	ofs(option->path_out);
			vcf_imputed->write(ofs, true);
		}
		else {
			ofstream	ofs(option->path_out, ios_base::app);
			vcf_imputed->write(ofs, false);		// headerを書かない
		}
		first = false;
		delete vcf_imputed;
	}
	
	delete vcf;
	delete sample_man;
	delete materials;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(1);
	}
	
	impute_VCF(option);
	delete option;
}