#ifndef __IMPUTEVCF
#define __IMPUTEVCF

#include <map>
#include <set>
#include "VCFHeteroHomo.h"
#include "VCFFillable.h"
#include "Pedigree.h"
#include "Map.h"

class Option;
class SampleManager;
class VCFFamily;
class VCFHeteroHomo;


//////////////////// Materials ////////////////////

class Materials {
	const std::string	path_map;
	const Map	*geno_map;
	const std::vector<const Map *>	chr_maps;
	
public:
	Materials(const std::string& path, const Map *m);
	~Materials();
	
	const std::string& get_path_map() const { return path_map; }
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const;
	double total_cM() const;
	void display_map_info() const;
	
public:
	static Materials *create(const Option *option);
};


//////////////////// process ////////////////////

VCFRecord *merge_progeny_records(std::vector<VCFFillable *>& vcfs,
									std::size_t i, const STRVEC& samples);
VCFSmall *impute_vcf_by_parents(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				const std::vector<const Family *>& families,
				const Map& geno_map, int num_threads);
VCFSmall *impute_vcf_by_parent(
				const VCFSmall *orig_vcf, VCFSmall *merged_vcf,
				const std::vector<const Family *>& families,
				const Map& geno_map,
				SampleManager *sample_man, int num_threads);
VCFSmall *impute_one_parent_vcf(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const std::vector<const Family *>& families,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads);
VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const std::vector<const Family *>& families,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads);
VCFSmall *impute_iolated_samples(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				SampleManager *sample_man, const STRVEC& samples,
				const Map& gmap, int num_threads);

void display_chromosome_info(const VCFSmall *orig_vcf);
VCFSmall *impute_vcf_chr(const VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option);
void print_info(const Option *option);
void impute_VCF(const Option *option);
#endif
