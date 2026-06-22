#ifndef __PARENTPROGENYIMPUTEDFAMILY
#define __PARENTPROGENYIMPUTEDFAMILY

#include <vector>
#include <string>

class VCFSmall;
class VCFGenoBase;
class VCFGeno;
class VCFImputable;
class VCFFamilyRecord;
class KnownFamily;
class Map;
class OptionSmall;


//////////////////// ParentProgenyImputedFamily ////////////////////

namespace ParentProgenyImputedFamily {
	bool is_small(const std::vector<std::vector<int>>& ref_haps,
					std::size_t L, std::size_t P, const OptionSmall& op);
	std::size_t compute_upper_num_progenies(
					const std::vector<std::vector<int>>& ref_haps,
					std::size_t L, std::size_t P, const OptionSmall& op);
	const KnownFamily *reduce_progenies(std::size_t p, const KnownFamily *f);
	VCFImputable *create_family_vcf(
							const KnownFamily *family,
							const std::vector<VCFFamilyRecord *>& records,
							std::size_t num_families,
							const std::vector<std::vector<int>>& ref_haps,
							bool should_impute_mat, const VCFSmall *orig_vcf,
							const OptionSmall& op);
	VCFGeno *impute(const VCFSmall *orig_vcf,
					const VCFGenoBase *imputed_vcf,
					const std::vector<const KnownFamily *>& families,
					const std::vector<std::string>& non_imputed_parents,
					const std::vector<std::vector<int>>& ref_haps,
					const OptionSmall& op);
};
#endif
