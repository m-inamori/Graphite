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
	const Map	*geno_map;
	const std::vector<const Map *>	chr_maps;
	
public:
	Materials(const Map *m);
	~Materials();
	
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const;
	double total_cM() const;
	void display_map_info() const;
	
public:
	static Materials *create(const Option *option);
};


//////////////////// process ////////////////////

typedef std::pair<std::string,std::string>	Parents;
typedef std::map<Parents, std::vector<VCFHeteroHomoRecord *>>	HeHoRecords;
typedef std::map<Parents, std::vector<VCFImpFamilyRecord *>>	ImpRecords;
typedef std::map<std::pair<Parents,bool>, VCFHeteroHomo *>	HeteroParentVCFs;

std::pair<HeHoRecords, ImpRecords> classify_records(VCFSmall *vcf,
						const std::vector<const Family *>& families, double p);
std::vector<std::vector<VCFImpFamilyRecord *>>
							sort_records(const HeHoRecords& heho_records,
											const ImpRecords& other_records);
void modify_00x11(const HeHoRecords& heho_records, ImpRecords& other_records);

std::pair<std::map<Parents, std::vector<VCFHeteroHomo *>>, HeHoRecords>
impute_hetero_homo_core(const HeHoRecords& records,
						const std::vector<const Family *>& families,
						const VCFSmall *vcf,
						const Map& geno_map, const Option *option);

void fill_in_thread(void *config);
void join_records(ImpRecords& records, HeHoRecords& unused_records);
std::vector<VCFFillable *> fill(
				std::map<Parents, std::vector<VCFHeteroHomo *>>& imputed_vcfs,
				ImpRecords& other_records, const Option *op);

HeteroParentVCFs extract_VCFs(const Materials *mat, const Option *option);
VCFFamily *impute_each(const Parents& parents, const Map& gmap,
										const HeteroParentVCFs& vcfs, int T);
VCFSmall *impute_vcf_chr(const VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option);
void impute_VCF(const Option *option);
#endif
