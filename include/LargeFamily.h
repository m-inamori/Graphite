#ifndef __LARGEFAMILY
#define __LARGEFAMILY

#include <vector>
#include "../include/VCF.h"

class VCFSmall;
class VCFFamily;
class VCFHeteroHomo;
class VCFFillable;
class VCFHeteroHomoRecord;
class VCFImpFamilyRecord;
class Family;
class Map;
class Option;
class TypeDeterminer;

namespace LargeFamily {
	typedef std::pair<std::string,std::string>	Parents;
	typedef std::map<Parents, std::vector<VCFHeteroHomoRecord *>>	HeHoRecords;
	typedef std::map<Parents, std::vector<VCFImpFamilyRecord *>>	ImpRecords;
	
	struct ConfigThreadClassify {
		std::size_t	first;
		std::size_t	num_threads;
		const VCFFamily	*vcf;
		const TypeDeterminer	*td;
		std::vector<VCFHeteroHomoRecord *>&	heho_records;
		std::vector<VCFImpFamilyRecord *>&	other_records;
		
		ConfigThreadClassify(std::size_t i, std::size_t	n,
							 const VCFFamily *vcf_, const TypeDeterminer *td_,
							 std::vector<VCFHeteroHomoRecord *>& heho_rs,
							 std::vector<VCFImpFamilyRecord *>& other_rs) :
					first(i), num_threads(n), vcf(vcf_), td(td_),
					heho_records(heho_rs), other_records(other_rs) { }
	};
	
	// HeteroHomoとそれ以外を分ける
	void classify_record(std::size_t i, const VCFFamily *vcf,
						 const TypeDeterminer *td,
						 std::vector<VCFHeteroHomoRecord *>& heho_records,
						 std::vector<VCFImpFamilyRecord *>& other_records);
	
	void classify_records_in_thread(void *config);
	
	// HeteroHomoだけ別にする
	// このあとHeteroHomoだけ補完するから
	// その他はVCFFillableにした後補完する
	std::pair<std::vector<VCFHeteroHomoRecord *>,
			  std::vector<VCFImpFamilyRecord *>>
	classify_records(const VCFFamily *vcf, const Option *option);
	
	struct ConfigThreadCreate {
		std::size_t	first;
		std::size_t	num_threads;
		const VCFSmall	*orig_vcf;
		const std::vector<const Family *>&	families;
		std::vector<VCFFamily *>&	vcfs_family;
		
		ConfigThreadCreate(std::size_t i, std::size_t n,
						   const VCFSmall *vcf,
						   const std::vector<const Family *>& fams,
						   std::vector<VCFFamily *>& vcfs) :
										first(i), num_threads(n), orig_vcf(vcf),
										families(fams), vcfs_family(vcfs) { }
	};
	
	void create_vcf_in_thread(void *config);
	std::vector<VCFFamily *> create_family_vcfs(
							const VCFSmall *orig_vcf,
							const std::vector<const Family *>& large_families,
							int num_threads);
	
	std::vector<std::pair<std::string,
						  std::vector<std::pair<std::size_t, std::size_t>>>>
	collect_same_parent_families(const std::vector<const Family *>& families);
	std::vector<std::vector<VCFImpFamilyRecord *>>
	collect_same_position_records(
				const std::vector<std::tuple<std::vector<VCFHeteroHomoRecord *>,
											 std::vector<VCFImpFamilyRecord *>,
											 int>>& rs);
	void modify_00x11_each(
				const std::vector<std::tuple<std::vector<VCFHeteroHomoRecord *>,
											 std::vector<VCFImpFamilyRecord *>,
											 int>>& rs);
	void modify_00x11(
				std::vector<std::vector<VCFHeteroHomoRecord *>>& heho_recordss,
				std::vector<std::vector<VCFImpFamilyRecord *>>& other_recordss,
				const std::vector<const Family *>& families);
	std::pair<std::vector<std::vector<VCFHeteroHomoRecord *>>,
			  std::vector<std::vector<VCFImpFamilyRecord *>>>
	divide_vcf_into_record_types(const std::vector<VCFFamily *>& family_vcfs,
														const Option *option);
	std::vector<VCFFillable *> fill_vcf(
			const std::map<std::string, std::vector<VCFHeteroHomo *>>& dic_vcfs,
			const std::vector<std::vector<VCFImpFamilyRecord *>>& other_recordss,
			const std::vector<const Family *>& families);
	void compress_records(std::vector<VCFImpFamilyRecord *>& others);
	VCFSmall *correct_large_family_VCFs(const VCFSmall *orig_vcf,
							const std::vector<const Family *>& large_families,
							const Map& geno_map, const Option *option);
}

#endif
