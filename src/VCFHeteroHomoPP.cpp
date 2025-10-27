#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFFillable.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFHeteroHomoPP ////////////////////

VCFHeteroHomoPP::VCFHeteroHomoPP(const STRVEC& s,
									const vector<VCFFillableRecord *>& rs,
									const Map& m, const VCFSmall *vcf) :
							VCFFamilyBase(s, vcf),
							VCFMeasurable(m), records(rs) { }

VCFHeteroHomoPP::~VCFHeteroHomoPP() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

bool VCFHeteroHomoPP::is_mat_hetero() const {
	return records.front()->is_hetero(0);
}

string VCFHeteroHomoPP::make_seq(size_t i) const {
	stringstream	ss;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		VCFFillableRecord	*record = *p;
		if(record->is_NA(i+2)) {
			ss << 'N';
			continue;
		}
		
		const int	a_mat1 = record->get_allele(0, 0);
		const int	a_pat1 = record->get_allele(1, 0);
		const int	a_pat2 = record->get_allele(1, 1);
		const int	prog_gt = record->unphased(i+2);
		if(a_pat1 == a_pat2) {	// mat hetero
			const int	diff = prog_gt - a_pat1;
			if(diff == -1 || diff == 2)
				ss << 'N';
			else
				ss << (a_mat1 == diff ? '0' : '1');
		}
		else {				// pat hetero
			const int	diff = prog_gt - a_mat1;
			if(diff == -1 || diff == 2)
				ss << 'N';
			else
				ss << (a_pat1 == diff ? '0' : '1');
		}
	}
	return ss.str();
}

string VCFHeteroHomoPP::impute_sample_seq(size_t i,
								const vector<double>& cMs, double min_c) {
	const string	seq = this->make_seq(i);
	if(Imputer::is_all_same_without_N(seq))
		return Imputer::create_same_color_string(seq, '0');
	
	const string	hidden_seq = Imputer::impute(seq, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

int VCFHeteroHomoPP::update_each(size_t i, size_t j, char c) {
	VCFFillableRecord	*record = this->records[i];
	const int	k = c == '0' ? 0 : 1;
	if(record->is_mat_hetero()) {
		const int	a1 = record->get_allele(0, k);
		const int	a2 = record->get_allele(1, 1);
		return Genotype::from_alleles(a1, a2);
	}
	else {
		const int	a1 = record->get_allele(0, 0);
		const int	a2 = record->get_allele(1, k);
		return Genotype::from_alleles(a1, a2);
	}
}

void VCFHeteroHomoPP::update(size_t i, const STRVEC& seqs) {
	for(size_t j = 2; j < this->get_samples().size(); ++j) {
		const char	c = seqs[j-2].c_str()[i];
		this->records[i]->set_geno(j, this->update_each(i, j, c));
	}
}

void VCFHeteroHomoPP::impute() {
	if(this->size() == 0)
		return;
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 0; i < this->num_progenies(); ++i) {
		imputed_seqs.push_back(impute_sample_seq(i, cMs, 1.0));
	}
	
	for(size_t k = 0; k < this->size(); ++k) {
		this->update(k, imputed_seqs);
	}
}

void VCFHeteroHomoPP::fill() {
	const Groups	*groups = Groups::create(records);
	const auto	record_sets = groups->create_record_sets();
	for(auto p = record_sets.begin(); p != record_sets.end(); ++p) {
		impute_core(*p);
	}
	delete groups;
	Common::delete_all(record_sets);
}

void VCFHeteroHomoPP::impute_core(const RecordSet *record_set) {
	auto	record = record_set->record;
	if(record == NULL)
		return;
	
	for(size_t i = 2; i < record->num_samples(); ++i) {
		const int	mat_from = record_set->determine_mat_from(i);
		const int	pat_from = record_set->determine_pat_from(i);
		const int	a1 = record->get_allele(0, mat_from - 1);
		const int	a2 = record->get_allele(1, pat_from - 1);
		record->set_geno(i, Genotype::from_alleles(a1, a2));
	}
}

pair<ParentComb, FillType> VCFHeteroHomoPP::classify_record(
													VCFFamilyRecord *record) {
	const int	i = record->is_mat_hetero() ? 0 : 1;
	const int	j = record->is_pat_hetero() ? 0 : 1;
	if(i == 0 && j == 0) {
		return make_pair(ParentComb::P01x01, FillType::IMPUTABLE);
	}
	else if(i == 1 && j == 0) {
		if(record->is_00(0))
			return make_pair(ParentComb::P00x01, FillType::PAT);
		else
			return make_pair(ParentComb::P01x11, FillType::PAT);
	}
	else if(i == 0 && j == 1) {
		if(record->is_00(1))
			return make_pair(ParentComb::P00x01, FillType::MAT);
		else
			return make_pair(ParentComb::P01x11, FillType::MAT);
	}
	else {
		if(record->is_00(0) && record->is_00(1))
			return make_pair(ParentComb::P00x00, FillType::FILLED);
		else if(record->is_11(0) && record->is_11(1))
			return make_pair(ParentComb::P11x11, FillType::FILLED);
		else
			return make_pair(ParentComb::P00x11, FillType::FILLED);
	}
}

array<vector<VCFFillableRecord *>, 4> VCFHeteroHomoPP::classify_records(
									const STRVEC& samples,
									const vector<VCFFamilyRecord *>& records,
									const VCFSmall *ref_vcf) {
	const auto	cols = ref_vcf->extract_columns(samples);
	array<vector<VCFFillableRecord *>, 4>	rss;
	for(size_t index = 0; index < records.size(); ++index) {
		VCFFamilyRecord	*record = records[index];
		const auto	pair = classify_record(record);
		const auto	*ref_record = ref_vcf->get_record(index);
		const auto	probs = ref_record->parse_PL(record->get_geno(), cols);
		auto	*new_record = new VCFFillableRecord(record->get_pos(),
												record->get_geno(), index,
												pair.second, pair.first, probs);
		const size_t	type_index = static_cast<size_t>(pair.second);
		rss[type_index].push_back(new_record);
	}
	return rss;
}

VCFFillable *VCFHeteroHomoPP::merge_vcf(const VCFHeteroHomoPP *mat_vcf,
					const VCFHeteroHomoPP *pat_vcf,
					const vector<VCFFillableRecord *>& homohomo_records,
					const vector<VCFFillableRecord *>& heterohetero_records) {
	for(auto p = homohomo_records.begin(); p != homohomo_records.end(); ++p) {
		VCFFillableRecord	*record = *p;
		const int	a1 = record->get_allele(0, 0);
		const int	a2 = record->get_allele(1, 0);
		const int	gt = Genotype::from_alleles(a1, a2);
		for(size_t i = 2; i != mat_vcf->get_samples().size(); ++i)
			record->set_geno(i, gt);
	}
	
	const vector<VCFFillableRecord *>&	mat_records = mat_vcf->get_records();
	const vector<VCFFillableRecord *>&	pat_records = pat_vcf->get_records();
	vector<VCFFillableRecord *>	records(mat_records.begin(), mat_records.end());
	records.insert(records.begin(), pat_records.begin(), pat_records.end());
	records.insert(records.begin(), homohomo_records.begin(),
										homohomo_records.end());
	records.insert(records.begin(), heterohetero_records.begin(),
										heterohetero_records.end());
	std::sort(records.begin(), records.end(), 
				[](const VCFFillableRecord *lh, const VCFFillableRecord *rh)
				{ return lh->get_pos() < rh->get_pos(); });
	return new VCFFillable(mat_vcf->get_samples(), records,
												mat_vcf->get_ref_vcf());
}

VCFFillableRecord *VCFHeteroHomoPP::merge_record(const VCFRecord *record1,
												 const VCFRecord *record2,
												 const STRVEC& samples, int i,
												 const TypeDeterminer *td) {
	const int	pos = record1->pos();
	vector<int>	geno;
	const auto&	v1 = record1->get_v();
	for(size_t c = 9; c != v1.size(); ++c)
		geno.push_back(Genotype::all_gt_to_int(v1[c]));
	const auto&	v2 = record2->get_v();
	for(size_t c = 9; c != v2.size(); ++c)
		geno.push_back(Genotype::all_gt_to_int(v2[c]));
	
	auto	*record = new VCFFamilyRecord(pos, geno);
	const auto	pair1 = VCFHeteroHomoPP::classify_record(record);
	const ParentComb	pc = pair1.first;
	const FillType		type = pair1.second;
	vector<VCFRecord::Probs>	probs;
	for(auto p = geno.begin(); p != geno.end(); ++p)
		probs.push_back(VCFRecord::decide_PL_by_genotype(*p));
	return new VCFFillableRecord(pos, geno, i, type, pc, probs);
}

VCFFillableRecord *VCFHeteroHomoPP::fill_NA(VCFRecord *record1,
											const STRVEC& samples, int i) {
	const int	pos = record1->pos();
	const size_t	NA_len = samples.size() - record1->num_samples();
	vector<int>	geno;
	const auto&	v = record1->get_v();
	for(size_t c = 9; c != v.size(); ++c)
		geno.push_back(Genotype::all_gt_to_int(v[c]));
	for(size_t i = 0; i < NA_len; ++i)
		geno.push_back(Genotype::NA);
	
	vector<VCFRecord::Probs>	probs;
	for(auto p = geno.begin(); p != geno.end(); ++p)
		probs.push_back(VCFRecord::decide_PL_by_genotype(*p));
	return new VCFFillableRecord(pos, geno, i,
									FillType::UNABLE, ParentComb::PNA, probs);
}

VCFHeteroHomoPP *VCFHeteroHomoPP::merge(const VCFSmall *vcf_parents,
										const VCFSmall *vcf_progenies,
										const VCFSmall *orig_vcf,
										const STRVEC& samples,
										const Map& m, const Option *option) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const TypeDeterminer	*td = CR->get_TypeDeterminer(samples.size()-2,
																option->ratio);
	const auto	header = vcf_parents->trim_header(samples);
	vector<VCFFillableRecord *>	records;
	size_t	j = 0;
	for(size_t i = 0; i < vcf_parents->size(); ++i) {
		auto	*record1 = vcf_parents->get_record(i);
		VCFFillableRecord	*record;
		if(j == vcf_progenies->size()) {
			record = VCFHeteroHomoPP::fill_NA(record1, samples, i);
		}
		else {
			auto	*prog_record = vcf_progenies->get_record(j);
			if(record1->pos() == prog_record->pos()) {
				record = VCFHeteroHomoPP::merge_record(record1,
															prog_record,
															samples, i, td);
				++j;
			}
			else {
				record = VCFHeteroHomoPP::fill_NA(record1, samples, i);
			}
		}
		records.push_back(record);
	}
	return new VCFHeteroHomoPP(samples, records, m, orig_vcf);
}
