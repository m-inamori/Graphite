#ifndef __VCFFAMILY
#define __VCFFAMILY

#include "VCF.h"


//////////////////// VCFFamilyRecord ////////////////////

class VCFFamilyRecord : public VCFRecord {
public:
	VCFFamilyRecord(const STRVEC& v, const STRVEC& samples) :
											VCFRecord(v, samples) { }
	~VCFFamilyRecord() { }
	
	VCFFamilyRecord *copy() const;
	const std::string& mat_gt() const { return this->v[9]; }
	const std::string& pat_gt() const { return this->v[10]; }
	int mat_int_gt() const { return this->get_int_gt(0); }
	int pat_int_gt() const { return this->get_int_gt(1); }
	std::vector<int> progeny_gts() const;
	std::vector<int> get_progeny_int_gts() const;
	bool is_mat_homo() const { return is_homo(9); }
	bool is_pat_homo() const { return is_homo(10); }
	std::size_t num_progenies() const { return samples.size() - 2U; }
	std::tuple<int,int,int> count_gts() const;
	
	void set_mat_GT(const std::string& gt) { set_GT(9, gt); }
	void set_pat_GT(const std::string& gt) { set_GT(10, gt); }
	void set_mat_int_GT(int gt) { set_int_GT(9, gt); }
	void set_pat_int_GT(int gt) { set_int_GT(10, gt); }
	
	void set(const STRVEC& new_v);
};


//////////////////// VCFFamily ////////////////////

class VCFFamily : public VCFSmall {
protected:
	std::vector<VCFFamilyRecord *>	family_records;
	
public:
	VCFFamily(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFamilyRecord *> rs);
	~VCFFamily() { }
	
	void set_records(const std::vector<VCFFamilyRecord *>& rs);
	void set_records_base(const std::vector<VCFFamilyRecord *>& rs);
	
	const std::string& mat() const { return samples[0]; }
	const std::string& pat() const { return samples[1]; }
	VCFFamilyRecord *get_record(std::size_t i) const {
		return family_records[i];
	}
	
	bool is_all_hetero(bool is_mat) const;
	bool is_all_homo(bool is_mat) const;
	
	void update_genotypes(const std::vector<STRVEC>& GTs);
	
public:
	static VCFFamily *merge(const VCFFamily *vcf1, const VCFFamily *vcf2);
	// join the VCFs divided into chrmosome
	static VCFFamily *join(const std::vector<VCFFamily *>& vcfs);
};
#endif
