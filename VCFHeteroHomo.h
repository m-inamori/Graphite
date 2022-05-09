#ifndef __VCFHETEROHOMO
#define __VCFHETEROHOMO

#include "VCFFamily.h"

class Map;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFHeteroHomoRecord ////////////////////

class VCFHeteroHomoRecord : public VCFFamilyRecord {
	enum class SEGTYPE { HomoHetero, HeteroHetero, HeteroHomo, None };
	
public:
	VCFHeteroHomoRecord(const STRVEC& v, const STRVEC& s) :
											VCFFamilyRecord(v, s) { }
	
	VCFHeteroHomoRecord *copy() const;
	int mat_int_gt() const { return get_int_gt(0); }
	int pat_int_gt() const { return get_int_gt(1); }
	std::vector<int> get_int_gts() const;
	
	SEGTYPE segregation_type() const;
	bool is_hetero_and_homo(bool is_mat) const;
	std::vector<int> genotypes_from_hetero_parent(bool is_mat_hetero) const;
	bool is_valid_segregation(bool is_mat, double cM) const;
	
private:
	std::vector<std::vector<double>> make_probability_table() const;
	std::vector<int> count_numbers() const;
	std::vector<double> log_likelihood(
						const std::vector<std::vector<double>>& pss,
						const std::vector<int>& ns) const;
	bool is_Mendelian_segregation() const;
	
	int genotype_from_hetero_parent(int i, int homo_gt) const;
	
/*
public:
	static const VCFFamilyRecord *subset(const VCFRecord *record,
											const std::vector<int>& columns);
*/
};


//////////////////// VCFHeteroHomo ////////////////////

class VCFHeteroHomo : public VCFFamily {
public:
	typedef std::pair<std::string,std::string>	Parents;
	
protected:
	std::vector<VCFHeteroHomoRecord *>	records;
	const Map&	genetic_map;
	
public:
	VCFHeteroHomo(const std::vector<STRVEC>& h, const STRVEC& s,
						std::vector<VCFHeteroHomoRecord *> rs, const Map& m);
	~VCFHeteroHomo();
	
	VCFHeteroHomoRecord *get_record(std::size_t i) const { return records[i]; }
	const Map& get_map() const { return genetic_map; }
	double cM(std::size_t i) const;
	std::pair<int,bool> distance(const std::vector<int>& gts1,
							const std::vector<int>& gts2, int max_dist) const;
	std::vector<VCFHeteroHomo *> divide_into_chromosomes() const;
	
	void update_genotypes(const std::vector<STRVEC>& GT_table);
	
private:
	std::vector<VCFFamilyRecord *> to_VCFFamilyRecord(
										std::vector<VCFHeteroHomoRecord *>& rs);
	
public:
	static std::vector<int> select_columns(const Parents& parents,
											const VCFOriginal *vcf,
											const PedigreeTable& pedigree);
	static std::map<std::pair<Parents,bool>, VCFHeteroHomo *> create_vcfs(
										VCFOriginal *orig_vcf,
										const std::vector<Parents>& families,
										const PedigreeTable& pedigree,
										const Map& geno_map,
										bool debug);
};
#endif
