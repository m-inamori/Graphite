#ifndef __SELFFAMILY
#define __SELFFAMILY

#include <vector>
#include <set>
#include <string>

class VCFSmall;
class VCFGeno;
class VCFSelfImputable;
class Family;
class KnownFamily;
class OptionSmall;


//////////////////// SelfFamily ////////////////////

namespace SelfFamily {
	GenoRecord *create_record(const GenoRecord *record,
								std::vector<std::size_t>& indices);
	std::vector<size_t> collect_progeny_indices(
								const Family *family,
								const std::set<std::string>& imputed_samples);
	VCFSelfImputable *create_family_vcf(const VCFSmall *orig_vcf,
								const std::vector<GenoRecord *>& records,
								const std::vector<std::vector<int>>& ref_haps,
								const KnownFamily *family,
								const std::set<std::string>& imputed_samples,
								const OptionSmall& op);
	VCFGeno *impute(const VCFSmall *orig_vcf,
						 const VCFGeno *imputed_vcf,
						 const std::vector<std::vector<int>>& ref_haps,
						 const std::vector<const KnownFamily *>& families,
						 const std::vector<std::string>& imputed_samples,
						 const OptionSmall& op);
};
#endif
