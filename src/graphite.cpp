#include <iostream>
#include <cassert>
#include "../include/graphite.h"
#include "../include/SampleManager.h"
#include "../include/VCFOriginal.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/VCFJunkRecord.h"
#include "../include/VCFHomoHomo.h"
#include "../include/VCFHeteroHeteroLite.h"
#include "../include/VCFFillable.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFIsolated.h"
#include "../include/option.h"
#include "../include/common.h"

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
		const auto	t = VCFHeteroHomo::make_VCFHeteroHomo(heho_records,
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
	
	const auto	w = VCFHeteroHomo::impute_hetero_homo_all(vcfs, option);
	
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

VCFSmall *fill_and_merge_vcf(
						map<Parents, vector<VCFHeteroHomo *>>& imputed_vcfs,
						ImpRecords& other_records,
						const STRVEC& samples, const Option *option) {
	// Familyごとに残りのRecordをphaseする
	const auto	filled_vcfs = VCFFillable::fill_all(imputed_vcfs,
														other_records,
														option->num_threads);
	// FamilyごとのVCFを全て統合する
	// この部分はVCFFillableでなく、切り離して別の名前空間で行いたい
	VCFSmall	*vcf_integrated = VCFFillable::merge(filled_vcfs, samples);
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
				const vector<const Family *>& families,
				const Map& geno_map, int num_threads) {
	vector<VCFFillable *>	vcfs = VCFHeteroHomoPP::impute_vcfs(
														orig_vcf, merged_vcf,
														families, geno_map,
														num_threads);
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
				const vector<const Family *>& families, const Map& geno_map,
				const SampleManager *sample_man, int num_threads) {
	vector<pair<const Family *, bool>>	parent_imputed_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		const bool	is_mat_phased = sample_man->is_imputed(family->get_mat());
		parent_imputed_families.push_back(make_pair(family, is_mat_phased));
	}
	
	const vector<VCFFamily *>	vcfs = VCFOneParentPhased::impute_all_by_parent(
														orig_vcf, merged_vcf,
														parent_imputed_families,
														geno_map, num_threads);
	
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

VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const vector<const Family *>& families,
								  SampleManager *sample_man, int num_threads) {
	vector<pair<const Family *, vector<size_t>>>	progeny_imputed_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		vector<size_t>	ppi;	// phased progeny index
		const vector<string>&	samples = family->get_samples();
		for(size_t i = 0; i < samples.size(); ++i) {
			if(sample_man->is_imputed(samples[i]))
				ppi.push_back(i);
		}
		progeny_imputed_families.push_back(make_pair(family, ppi));
	}
	const auto	vcfs = VCFProgenyPhased::impute_all_by_progeny(orig_vcf,
													merged_vcf,
													progeny_imputed_families,
													num_threads);
	
	vector<VCFSmall *>	vcfs2(vcfs.begin(), vcfs.end());	// convert class
	VCFSmall	*new_vcf = VCFSmall::join(vcfs2, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}

VCFSmall *impute_iolated_samples(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				SampleManager *sample_man, const STRVEC& samples,
				const Map& gmap, int num_threads) {
	const STRVEC	references = sample_man->get_large_parents();
	// あとでマルチプロセス化するためにphasingすべきsample分割する
	auto	vcfs = VCFIsolated::create(orig_vcf, merged_vcf,
										samples, references, gmap, num_threads);
	const auto	new_vcfs = VCFIsolated::impute_all(vcfs, num_threads);
	VCFSmall	*new_vcf = VCFSmall::join(new_vcfs, orig_vcf->get_samples());
	Common::delete_all(new_vcfs);
	return new_vcf;
}

VCFSmall *impute_vcf_chr(const VCFSmall *orig_vcf, SampleManager *sample_man,
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
												merged_vcf, families1,
												geno_map, option->num_threads);
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
			auto	*new_imputed_vcf = impute_vcf_by_parent(orig_vcf,
											merged_vcf, families2, geno_map,
											sample_man, option->num_threads);
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
														families3, sample_man,
														option->num_threads);
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
	
	// 最後に孤立したサンプルを補完する
	const STRVEC	samples = sample_man->extract_isolated_samples();
	if(!samples.empty()) {
		VCFSmall	*new_imputed_vcf = impute_iolated_samples(
												orig_vcf, merged_vcf, sample_man,
												samples,
												geno_map, option->num_threads);
		vector<VCFSmall *>	vcfs{ merged_vcf, new_imputed_vcf };
		merged_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
		delete new_imputed_vcf;
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
