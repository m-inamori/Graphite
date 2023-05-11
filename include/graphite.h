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
	Materials(const Map *m) : geno_map(m),
						chr_maps(m->divide_into_chromosomes()) { }
	~Materials();
	
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const { return chr_maps[i]; }
	
public:
	static Materials *create(const Option *option);
};


//////////////////// process ////////////////////

typedef std::pair<std::string,std::string>	Parents;
typedef std::map<Parents, std::vector<VCFHeteroHomoRecord *>>	HeHoRecords;
typedef std::map<Parents, std::vector<VCFImpFamilyRecord *>>	ImpRecords;
typedef std::map<std::pair<Parents,bool>, VCFHeteroHomo *>	HeteroParentVCFs;
typedef std::pair<std::vector<VCFHeteroHomo *>,
					std::vector<VCFImpFamilyRecord *>>	Item;

std::pair<HeHoRecords, ImpRecords> classify_records(VCFSmall *vcf,
						const std::vector<const Family *>& families, double p);
std::vector<std::vector<VCFImpFamilyRecord *>>
							sort_records(const HeHoRecords& heho_records,
											const ImpRecords& other_records);
void modify_00x11(const HeHoRecords& heho_records, ImpRecords& other_records);

struct ConfigImpThread {
	const std::vector<std::vector<VCFHeteroHomo *>>&	vcfs_heho;
	const Option *option;
	const std::size_t	first;
	const int	num_threads;
	std::vector<ImpResult>&	imputed_vcfs;
	
	ConfigImpThread(const std::vector<std::vector<VCFHeteroHomo *>>& vcfs,
											const Option *op, int f, int n,
											std::vector<ImpResult>& results) :
									vcfs_heho(vcfs), option(op), first(f),
									num_threads(n), imputed_vcfs(results) { }
	
	std::size_t size() const { return imputed_vcfs.size(); }
};

std::vector<ImpResult> impute_hetero_homo_parellel(
				const std::map<std::string, std::vector<VCFHeteroHomo *>>& vcfs,
				const Option *option);

std::tuple<VCFHeteroHomo *, VCFHeteroHomo *, std::vector<VCFHeteroHomoRecord *>>
make_VCFHeteroHomo(const std::vector<VCFHeteroHomoRecord *>& records,
						const Family *family,
						const VCFSmall *vcf, const Map& geno_map);
std::pair<std::map<Parents, std::vector<VCFHeteroHomo *>>, HeHoRecords>
impute_hetero_homo_core(const HeHoRecords& records,
						const std::vector<const Family *>& families,
						const VCFSmall *vcf,
						const Map& geno_map, const Option *option);

struct ConfigFillThread {
	const std::vector<Item>&	items;
	const bool	all_out;
	const std::size_t	first;
	const int	num_threads;
	std::vector<VCFFillable *>&	filled_vcfs;
	
	ConfigFillThread(const std::vector<Item>& items_, bool ao, int f, int n,
										std::vector<VCFFillable *>& results) :
									items(items_), all_out(ao), first(f),
									num_threads(n), filled_vcfs(results) { }
	
	std::size_t size() const { return items.size(); }
};

void fill_in_thread(void *config);
std::vector<VCFFillable *> fill_parellel(
								std::vector<Item>& items, const Option *op);
void join_records(ImpRecords& records, HeHoRecords& unused_records);
std::vector<VCFFillable *> fill(
				std::map<Parents, std::vector<VCFHeteroHomo *>>& imputed_vcfs,
				ImpRecords& other_records, const Option *op);

HeteroParentVCFs extract_VCFs(const Materials *mat, const Option *option);
VCFFamily *impute_each(const Parents& parents, const Map& gmap,
										const HeteroParentVCFs& vcfs, int T);
VCFSmall *impute_vcf_chr(VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option);
void impute_VCF(const Option *option);
#endif
