#include <iostream>
#include <cassert>
#include "impute_VCF.h"
#include "VCFOriginal.h"
#include "VCFHeteroHomo.h"
#include "VCFJunkRecord.h"
#include "VCFHomoHomo.h"
#include "VCFHeteroHeteroLite.h"
#include "VCFFillable.h"
#include "option.h"
#include "common.h"

using namespace std;


//////////////////// Materials ////////////////////

Materials::~Materials() {
	delete pedigree;
	delete geno_map;
	Common::delete_all(chr_maps);
	Common::delete_all(families);
}

void Materials::select_families(set<pair<string,string>>& set_families) {
	vector<const Family *>	old_families = families;
	families.clear();
	for(auto p = old_families.begin(); p != old_families.end(); ++p) {
		if(set_families.find((*p)->parents()) != set_families.end())
			families.push_back(*p);
	}
}

vector<const Family *> Materials::make_families(
						const PedigreeTable *pedigree, const Option *option) {
	const vector<const Family *>	families = pedigree->extract_families();
	if(option->families.empty())
		return families;
	
	const set<int>	set_families(option->families.begin(),
								 option->families.end());
	vector<const Family *>	filtered_families;
	for(int i = 0; i < (int)families.size(); ++i) {
		if(set_families.find(i) != set_families.end())
			filtered_families.push_back(families[i]);
		else
			delete families[i];
	}
	return filtered_families;
}

Materials *Materials::create(const Option *option) {
	auto	*vcf = VCFOriginal::read(option->path_vcf);
	const auto	*pedigree = PedigreeTable::create(option->path_ped,
														vcf->get_samples());
	const auto	*geno_map = Map::read(option->path_map);
	auto	families = make_families(pedigree, option);
	delete vcf;
	return new Materials(pedigree, geno_map, families);
}


//////////////////// process ////////////////////

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

vector<pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>>
impute_hetero_homo_parellel(
			map<string, vector<VCFHeteroHomo *>>& vcfs, const Option *option) {
	vector<pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>>	v;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		vector<VCFHeteroHomo *>&	vcfs_heho = p->second;
		v.push_back(VCFHeteroHomo::impute_vcfs(vcfs_heho, option));
	}
	return v;
}

pair<map<Parents, vector<VCFHeteroHomo *>>, HeHoRecords> impute_hetero_homo(
								const HeHoRecords& records,
								const vector<const Family *>& families,
								const VCFSmall *vcf,
								const Map& geno_map, const Option *option) {
	// あとでFamilyを回復するために必要
	map<Parents, const Family *>	parents_to_family;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		parents_to_family.insert(make_pair(family->parents(), family));
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
		
		if(!option->all_out)
			continue;
		
		for(auto q = unused_records2.begin(); q != unused_records2.end(); ++q) {
			VCFHeteroHomoRecord	*record = *q;
			record->enable_modification();
			unused_records[record->parents()].push_back(record);
		}
	}
	return make_pair(imputed_vcfs, unused_records);
}

// HeteroHomoだけ別にする
// このあとHeteroHomoだけ補完するから
// その他はVCFFillableにした後補完する
pair<HeHoRecords, ImpRecords> classify_records(VCFSmall *vcf,
							const vector<const Family *>& families, double p) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	HeHoRecords	heho_records;
	ImpRecords	other_records;
	for(auto q = families.begin(); q != families.end(); ++q) {
		const Family	*family = *q;
		const Parents	parents = family->parents();
		const STRVEC&	samples = family->get_samples();
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

VCFFillable *fill_vcf(Item &v, bool all_out) {
	auto	vcfs = v.first;
	auto	records = v.second;
	return VCFFillable::fill(vcfs, records, all_out);
}

vector<VCFFillable *> fill(vector<Item>& items, bool all_out) {
	vector<VCFFillable *>	filled_vcfs;
	for(auto p = items.begin(); p != items.end(); ++p) {
		auto	*vcf = fill_vcf(*p, all_out);
		filled_vcfs.push_back(vcf);
	}
	return filled_vcfs;
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

VCFSmall *impute_vcf_chr(VCFSmall *vcf, const Map *geno_map,
				const vector<const Family *>& families, const Option *option) {
	cerr << "chr : " << vcf->get_records().front()->chrom() << endl;
	const auto	p = classify_records(vcf, families, option->ratio);
	const auto	heho_records = p.first;
	auto	other_records = p.second;
	modify_00x11(heho_records, other_records);
	
	// HeteroHomoとそれ以外を分けて、HeteroHomoを補完する
	// 使われなかったHeteroHomoはfillに回す
	const auto	p2 = impute_hetero_homo(heho_records, families,
											vcf, *geno_map, option);
	auto	imputed_vcfs = p2.first;
	auto	unused_records = p2.second;
	
	for(auto q = unused_records.begin(); q != unused_records.end(); ++q) {
		const auto&	parents = q->first;
		const auto&	rs = q->second;
		auto&	records = other_records[parents];
		records.insert(records.end(), rs.begin(), rs.end());
	}
	
	vector<Item>	items;
	for(auto q = imputed_vcfs.begin(); q != imputed_vcfs.end(); ++q) {
		const Parents&	parents = q->first;
		vector<VCFHeteroHomo *>&	vcfs = q->second;
		items.push_back(make_pair(vcfs, other_records[parents]));
	}
	vector<VCFFillable *>	vcfs_complemented = fill(items, option->all_out);
	for(auto q = items.begin(); q != items.end(); ++q) {
		for(auto r = q->first.begin(); r != q->first.end(); ++r)
			delete *r;
		for(auto r = q->second.begin(); r != q->second.end(); ++r)
			delete *r;
	}
	VCFSmall	*vcf_integrated = VCFFillable::merge(vcfs_complemented,
												vcf->get_samples(), option);
	for(auto p = vcfs_complemented.begin(); p != vcfs_complemented.end(); ++p)
		delete *p;
	return vcf_integrated;
}

void impute_VCF(const Option *option) {
	// chromosomeごとに処理する
	Materials	*materials = Materials::create(option);
	VCFHuge	*vcf = VCFHuge::read(option->path_vcf);
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
		const VCFSmall	*vcf_imputed = impute_vcf_chr(vcf_chrom, gmap,
											materials->get_families(), option);
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
