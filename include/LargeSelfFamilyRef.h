#ifndef __LARGEFAMILY
#define __LARGEFAMILY

#include <vector>
#include <set>

class VCFSmall;
class VCFGeno;
class VCFSelfFillable;
class VCFSelfParentImputed;
class VCFSelfProgenyImputed;
class KnownFamily;
class Map;
class Option;

namespace LargeSelfFamilyRef {
	VCFSelfParentImputed *create_parent_phased_vcf(const KnownFamily *family,
												   const VCFSmall *orig_vcf,
												   const VCFGeno *phased_vcf,
												   const Map& gmap,
												   const Option& op);
	VCFSelfProgenyImputed *create_progeny_phased_vcf(const KnownFamily *family,
													 const VCFSmall *orig_vcf,
													 const VCFGeno *ref_vcf,
													 const VCFGeno *phased_vcf,
													 const Map& gmap,
													 const Option& op);
	VCFGeno *extract_parents(const std::vector<VCFSelfFillable *>& vcfs);
	bool is_intersect(const std::vector<std::string>& v,
							const std::set<std::string>& s);
	VCFGeno *impute(const std::vector<const KnownFamily *>& families,
								const VCFSmall *orig_vcf,
								VCFGeno *merged_vcf,
								const VCFGeno *ref_vcf,
								const Map& geno_map, const Option& op);
}
#endif
