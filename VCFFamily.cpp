#include <algorithm>
#include <cassert>
#include "VCFFamily.h"

using namespace std;


//////////////////// VCFFamilyRecord ////////////////////

VCFFamilyRecord *VCFFamilyRecord::copy() const {
	return new VCFFamilyRecord(v, samples);
}

vector<int> VCFFamilyRecord::progeny_gts() const {
	vector<int>	w(this->num_progenies());
	for(size_t i = 2U; i < samples.size(); ++i)
		w[i-2] = this->get_int_gt(i);
	return w;
}

vector<int> VCFFamilyRecord::get_progeny_int_gts() const {
	vector<int>	w(this->num_progenies());
	for(size_t i = 2U; i < this->samples.size(); ++i) {
		assert(i < w.size() + 2);
		w[i-2] = get_int_gt(i);
	}
	return w;
}

tuple<int,int,int> VCFFamilyRecord::count_gts() const {
	int	counter[] = { 0, 0, 0 };
	const auto	gts = progeny_gts();
	for(auto p = gts.begin(); p != gts.end(); ++p) {
		if(0 <= *p && *p <= 2)
			counter[*p] += 1;
	}
	return make_tuple(counter[0], counter[1], counter[2]);
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

void VCFFamily::set_records(const vector<VCFFamilyRecord *>& rs) {
	family_records = rs;
	// ここに実装を書くと、深い継承の場合、最上位まで全部書かなけらばならない
	// 連鎖的にする
	set_records_base(rs);
}

// 親に書くと子クラスのvectorを全て並べることになるので子に書く
void VCFFamily::set_records_base(const vector<VCFFamilyRecord *>& rs) {
	records.insert(records.end(), rs.begin(), rs.end());
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

VCFFamily *VCFFamily::create(const VCFSmall *vcf, const STRVEC& samples) {
	const auto	columns = VCFFamily::select_columns(samples, vcf);
	if(columns.size() < 10)
		return NULL;
	
	const auto	header = vcf->create_header(samples);
	const vector<VCFRecord *>&	orig_records = vcf->get_records();
	vector<VCFFamilyRecord *>	records;
	for(auto p = orig_records.begin(); p != orig_records.end(); ++p)
		records.push_back(VCFFamily::subset(*p, samples, columns));
	VCFFamily	*new_vcf = new VCFFamily(header, samples, records);
	return new_vcf;
}

vector<int> VCFFamily::select_columns(const STRVEC& samples,
												const VCFSmall *vcf) {
	const STRVEC&	orig_samples = vcf->get_samples();
	map<string, int>	dic;
	int	c = 9;
	for(auto p = orig_samples.begin(); p != orig_samples.end(); ++p, ++c)
		dic.insert(make_pair(*p, c));
	
	vector<int>	columns;
	for(auto p = samples.begin(); p != samples.end(); ++p)
		columns.push_back(dic[*p]);
	return columns;
}

VCFFamilyRecord *VCFFamily::subset(VCFRecord *record, const STRVEC& samples,
												const vector<int>& columns) {
	const STRVEC&	orig_v = record->get_v();
	STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
	for(auto p = columns.begin(); p != columns.end(); ++p) {
		v.push_back(orig_v[*p]);
	}
	return new VCFFamilyRecord(v, samples);
}

VCFFamily *VCFFamily::merge(const VCFFamily *vcf1, const VCFFamily *vcf2) {
	// samplesをVFCとRecordで共通にするために
	// 空のVCFを作って、そのsamplesを使うようにする
	vector<VCFFamilyRecord *>	records;
	VCFFamily	*vcf = new VCFFamily(vcf1->get_header(),
											vcf1->get_samples(), records);
	vcf1->copy_chrs(vcf);
	
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
	
	vcf->set_records(records);
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
