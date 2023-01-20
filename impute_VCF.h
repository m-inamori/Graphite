#ifndef __IMPUTEVCF
#define __IMPUTEVCF

#include <map>
#include <set>
#include "VCFHeteroHomo.h"
#include "VCFFillable.h"
#include "Pedigree.h"
#include "Map.h"

class Option;
class VCFFamily;
class VCFHeteroHomo;


//////////////////// Materials ////////////////////

class Materials {
	const PedigreeTable	*pedigree;
	const Map	*geno_map;
	const std::vector<const Map *>	chr_maps;
	std::vector<const Family *>	families;
	
public:
	Materials(const PedigreeTable *p, const Map *m,
							const std::vector<const Family *>& f) :
						pedigree(p), geno_map(m),
						chr_maps(m->divide_into_chromosomes()), families(f) { }
	~Materials();
	
	const PedigreeTable& get_ped() const { return *pedigree; }
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const { return chr_maps[i]; }
	const std::vector<const Family *>& get_families() const { return families; }
	
	void select_families(
				std::set<std::pair<std::string,std::string>>& set_families);
	
public:
	static std::vector<const Family *> make_families(
						const PedigreeTable *pedigree, const Option *option);
	static Materials *create(const Option *option);
};


//////////////////// process ////////////////////

typedef std::pair<std::string,std::string>	Parents;
typedef std::map<Parents, std::vector<VCFHeteroHomoRecord *>>	HeHoRecords;
typedef std::map<Parents, std::vector<VCFImpFamilyRecord *>>	ImpRecords;
typedef std::map<std::pair<Parents,bool>, VCFHeteroHomo *>	HeteroParentVCFs;
typedef std::pair<std::vector<VCFHeteroHomo *>,
					std::vector<VCFImpFamilyRecord *>>	Item;

struct ConfigThread {
	const std::vector<Parents>&	families;
	const HeteroParentVCFs&	vcfs;
	const std::vector<const Map *>&	chr_maps;
	const std::size_t	first;
	const int	num_threads;
	std::vector<VCFFamily *>& imputed_vcfs;
	
	ConfigThread(const std::vector<Parents>& fams,
					const HeteroParentVCFs& v,
					const std::vector<const Map *>& m,
					int f, int n, std::vector<VCFFamily *>& results) :
								families(fams), vcfs(v), chr_maps(m), first(f),
								num_threads(n), imputed_vcfs(results) { }
};

std::vector<std::pair<std::vector<VCFHeteroHomo *>,
						std::vector<VCFHeteroHomoRecord *>>>
		impute_hetero_homo_parellel(
				std::map<std::string, std::vector<VCFHeteroHomo *>>& vcfs,
				const Option *option);

std::tuple<VCFHeteroHomo *, VCFHeteroHomo *, std::vector<VCFHeteroHomoRecord *>>
make_VCFHeteroHomo(const std::vector<VCFHeteroHomoRecord *>& records,
						const Family *family,
						const VCFSmall *vcf, const Map& geno_map);
std::vector<std::vector<VCFImpFamilyRecord *>>
							sort_records(const HeHoRecords& heho_records,
											const ImpRecords& other_records);
void modify_00x11(const HeHoRecords& heho_records, ImpRecords& other_records);
std::pair<std::map<Parents, std::vector<VCFHeteroHomo *>>, HeHoRecords>
impute_hetero_homo(const HeHoRecords& records,
					const std::vector<const Family *>& families,
					const VCFSmall *vcf,
					const Map& geno_map, const Option *option);
VCFFillable *fill_vcf(Item &v, bool all_out);
std::vector<VCFFillable *> fill(std::vector<Item>& items, bool all_out);

void select_families(Materials *materials, const HeteroParentVCFs& vcfs);
HeteroParentVCFs extract_VCFs(const Materials *mat, const Option *option);
VCFFamily *impute_each(const Parents& parents, const Map& gmap,
										const HeteroParentVCFs& vcfs, int T);
VCFSmall *impute_vcf_chr(VCFSmall *vcf,
							const std::vector<const Family *>& families,
							const Map *geno_map, const Option *option);
void impute_VCF(const Option *option);
#endif
