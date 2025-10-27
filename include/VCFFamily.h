#ifndef __VCFFAMILY
#define __VCFFAMILY

#include "GenoRecord.h"
#include "VCFGeno.h"


using Parents = std::pair<std::string, std::string>;


//////////////////// VCFFamilyRecord ////////////////////

class VCFFamilyRecord : public GenoRecord {
public:
	VCFFamilyRecord(const ll pos, const std::vector<int>& geno) :
											GenoRecord(pos, geno) { }
	~VCFFamilyRecord() { }
	
	VCFFamilyRecord *copy() const;
	int mat_gt() const { return geno[0]; }
	int pat_gt() const { return geno[1]; }
	int unphased_mat() const { return unphased(0); }
	int unphased_pat() const { return unphased(1); }
	bool is_mat_homo() const { return is_homo(0); }
	bool is_pat_homo() const { return is_homo(1); }
	bool is_mat_ref_homo() const { return is_ref_homo(0); }
	bool is_pat_ref_homo() const { return is_ref_homo(1); }
	bool is_mat_alt_homo() const { return is_alt_homo(0); }
	bool is_pat_alt_homo() const { return is_alt_homo(1); }
	bool is_mat_hetero() const { return is_hetero(0); }
	bool is_pat_hetero() const { return is_hetero(1); }
	bool is_mat_NA() const { return is_NA(0); }
	bool is_pat_NA() const { return is_NA(1); }
	
	int get_mat_allele(int j) const { return get_allele(0, j); }
	int get_pat_allele(int j) const { return get_allele(1, j); }
	
	std::vector<int> get_progeny_int_gts() const;
	std::vector<int> progeny_gts() const;
	std::size_t num_progenies() const { return geno.size() - 2; }
};


//////////////////// VCFFamilyBase ////////////////////

class VCFFamilyBase : public VCFGenoBase {
public:
	VCFFamilyBase(const STRVEC& s, const VCFSmall *vcf) :
											VCFGenoBase(s, vcf) { }
	virtual ~VCFFamilyBase() { };
	
	virtual VCFFamilyRecord *get_family_record(std::size_t i) const = 0;
	const std::string& mat() const { return samples[0]; }
	const std::string& pat() const { return samples[1]; }
	std::size_t num_progenies() const { return samples.size() - 2; }
	
	const Parents parents() const { return Parents(mat(), pat()); }
};


//////////////////// VCFFamily ////////////////////

class VCFFamily : public VCFFamilyBase {
protected:
	std::vector<VCFFamilyRecord *>	records;
	
public:
	VCFFamily(const STRVEC& s, const std::vector<VCFFamilyRecord *>& rs,
														const VCFSmall *vcf);
	~VCFFamily();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	const std::vector<VCFFamilyRecord *>& get_family_records() const {
		return records;
	}
	void set_records(const std::vector<VCFFamilyRecord *>& rs) { records = rs; }
	void clear_records() { records.clear(); }
	
public:
	static VCFFamily *create(const VCFSmall *vcf, const STRVEC& samples);
	static VCFFamilyRecord *subset(VCFRecord *record, const STRVEC& samples,
									const std::vector<std::size_t>& columns);
	// Create a new VCF by obtaining the Genotypes of the parents from vcf1
	// and the progenies from vcf2
	static VCFFamily *create_by_two_vcfs(const VCFGenoBase *vcf1,
										 const VCFSmall *vcf2,
										 const STRVEC& samples);
	static VCFFamilyRecord *subset(const VCFRecord *record,
									const std::vector<std::size_t>& columns);
	static VCFFamily *convert(const VCFFamilyBase *vcf);
};
#endif
