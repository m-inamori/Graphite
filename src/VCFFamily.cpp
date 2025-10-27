#include <algorithm>
#include <cassert>
#include "../include/VCFFamily.h"

using namespace std;


//////////////////// VCFFamilyRecord ////////////////////

VCFFamilyRecord *VCFFamilyRecord::copy() const {
	return new VCFFamilyRecord(pos, geno);
}

vector<int> VCFFamilyRecord::get_progeny_int_gts() const {
	vector<int>	w(this->num_progenies());
	for(size_t i = 0; i < this->num_progenies(); ++i)
		w[i] = this->unphased(i+2);
	return w;
}

vector<int> VCFFamilyRecord::progeny_gts() const {
	vector<int>	w(this->num_progenies());
	std::copy(geno.begin() + 2, geno.end(), w.begin());
	return w;
}


//////////////////// VCFFamily ////////////////////

VCFFamily::VCFFamily(const STRVEC& s, const vector<VCFFamilyRecord *>& rs,
														const VCFSmall *vcf) :
										VCFFamilyBase(s, vcf), records(rs) { }

VCFFamily::~VCFFamily() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

VCFFamily *VCFFamily::create(const VCFSmall *vcf, const STRVEC& samples) {
	const auto	columns = vcf->extract_columns(samples);
	const vector<VCFRecord *>&	orig_records = vcf->get_records();
	vector<VCFFamilyRecord *>	records;
	for(auto p = orig_records.begin(); p != orig_records.end(); ++p)
		records.push_back(VCFFamily::subset(*p, samples, columns));
	return new VCFFamily(samples, records, vcf);
}

VCFFamilyRecord *VCFFamily::subset(VCFRecord *record, const STRVEC& samples,
												const vector<size_t>& columns) {
	const STRVEC&	orig_v = record->get_v();
	const ll	pos = record->pos();
	vector<int>	geno(columns.size());
	for(size_t i = 0; i < columns.size(); ++i) {
		const size_t	c = columns[i];
		if(c != string::npos)
			geno[i] = Genotype::all_gt_to_int(orig_v[c]);
		else
			geno[i] = Genotype::NA;
	}
	return new VCFFamilyRecord(pos, geno);
}

VCFFamily *VCFFamily::create_by_two_vcfs(const VCFGenoBase *vcf1,
										 const VCFSmall *vcf2,
										 const STRVEC& samples) {
	// assume that samples are [mat, pat, prog1, prog2, ...]
	const vector<size_t>	columns1 = vcf1->extract_columns(samples);
	const vector<size_t>	columns2 = vcf2->extract_columns(samples);
	
	vector<VCFFamilyRecord *>	new_records;
	for(size_t i = 0; i < vcf1->size(); ++i) {
		const ll	pos = vcf1->get_pos(i);
		const auto&	geno1 = vcf1->get_record(i)->get_geno();
		const auto&	v2 = vcf2->get_record(i)->get_v();
		vector<int>	geno;
		for(size_t j = 0; j < samples.size(); ++j) {
			if(columns1[j] != string::npos)
				geno.push_back(geno1[columns1[j]]);
			else if(columns2[j] != string::npos)
				geno.push_back(Genotype::all_gt_to_int(v2[columns2[j]]));
			else
				geno.push_back(Genotype::NA);
		}
		auto	*new_record = new VCFFamilyRecord(pos, geno);
		new_records.push_back(new_record);
	}
	return new VCFFamily(samples, new_records, vcf1->get_ref_vcf());
}

VCFFamilyRecord *VCFFamily::subset(const VCFRecord *record,
									const vector<size_t>& columns) {
	const auto&	v = record->get_v();
	const ll	pos = record->pos();
	vector<int>	geno;
	for(auto p = columns.begin(); p != columns.end(); ++p) {
		if(*p != string::npos)
			geno.push_back(Genotype::all_gt_to_int(v[*p]));
		else
			geno.push_back(Genotype::NA);
	}
	return new VCFFamilyRecord(pos, geno);
}

VCFFamily *VCFFamily::convert(const VCFFamilyBase *vcf) {
	vector<VCFFamilyRecord *>	records(vcf->size());
	for(size_t i = 0; i < vcf->size(); ++i) {
		records[i] = vcf->get_family_record(i);
	}
	return new VCFFamily(vcf->get_samples(), records, vcf->get_ref_vcf());
}
