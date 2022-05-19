#include <algorithm>
#include <cassert>
#include "VCFFamily.h"

using namespace std;


//////////////////// VCFFamilyRecord ////////////////////

VCFFamilyRecord *VCFFamilyRecord::copy() const {
	return new VCFFamilyRecord(v, samples);
}

vector<int> VCFFamilyRecord::progeny_gts() const {
	vector<int>	w(this->samples.size() - 2);
	for(size_t i = 2U; i < samples.size(); ++i)
		w[i-2] = this->get_int_gt(i);
	return w;
}

bool VCFFamilyRecord::is_homo(size_t i) const {
	const string&	s = this->v[i+9];
	return s.c_str()[0] == s.c_str()[2];
}

void VCFFamilyRecord::set(const STRVEC& new_v) {
	v = new_v;
	VCFRecord::set(v);
}


//////////////////// VCFFamily ////////////////////

VCFFamily::VCFFamily(const vector<STRVEC>& h, const STRVEC& s,
									vector<VCFFamilyRecord *> rs) :
					VCFSmall(h, s, vector<VCFRecord *>(rs.begin(), rs.end())),
					family_records(rs) {
}

bool VCFFamily::is_all_hetero(bool is_mat) const {
	for(auto p = family_records.begin(); p != family_records.end(); ++p) {
		if(is_mat) {
			if((*p)->mat_int_gt() != 1)
				return false;
		}
		else {
			if((*p)->pat_int_gt() != 1)
				return false;
		}
	}
	return true;
}

bool VCFFamily::is_all_homo(bool is_mat) const {
	for(auto p = family_records.begin(); p != family_records.end(); ++p) {
		if(is_mat) {
			if((*p)->mat_int_gt() == 1)
				return false;
		}
		else {
			if((*p)->pat_int_gt() == 1)
				return false;
		}
	}
	return true;
}

void VCFFamily::update_genotypes(const std::vector<STRVEC>& GTs) {
	for(size_t i = 0U; i < records.size(); ++i)
		records[i]->set_GTs(GTs[i]);
	VCFSmall::update_genotypes(GTs);
}

VCFFamily *VCFFamily::merge(const VCFFamily *vcf1, const VCFFamily *vcf2) {
	vector<VCFFamilyRecord *>	records;
	size_t	k = 0U;
	size_t	l = 0U;
	while(k < vcf1->size() && l < vcf2->size()) {
		VCFFamilyRecord	*record1 = vcf1->get_record(k);
		VCFFamilyRecord	*record2 = vcf2->get_record(l);
		POSITION	pos1 = vcf1->record_position(*record1);
		POSITION	pos2 = vcf2->record_position(*record2);
		assert(pos1 != pos2);
		if(pos1 < pos2) {
			records.push_back(record1->copy());
			k += 1;
		}
		else {
			records.push_back(record2->copy());
			l += 1;
		}
	}
	
	for(size_t i = k; i < vcf1->size(); ++i)
		records.push_back(vcf1->get_record(i)->copy());
	for(size_t i = l; i < vcf2->size(); ++i)
		records.push_back(vcf2->get_record(i)->copy());
	
	VCFFamily	*vcf = new VCFFamily(vcf1->get_header(),
										vcf1->get_samples(), records);
	vcf1->copy_chrs(vcf);
	return vcf;
}

VCFFamily *VCFFamily::join(const vector<VCFFamily *>& vcfs) {
	vector<VCFFamilyRecord *>	records;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFamily	*vcf = *p;
		for(size_t i = 0U; i < vcf->size(); ++i) {
			VCFFamilyRecord	*record = vcf->get_record(i);
			records.push_back(record->copy());	// for memory leak
		}
	}
	VCFFamily	*vcf = vcfs.front();
	return new VCFFamily(vcf->get_header(), vcf->get_samples(), records);
}
