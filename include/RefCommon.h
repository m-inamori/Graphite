#ifndef __REFCOMMON
#define __REFCOMMON

#include <vector>
#include "VCF.h"

class VCFFamily;
class GenoRecord;
class VCFGenoBase;
class VCFGeno;
class VCFFamilyRecord;
class KnownFamily;
class Map;


//////////////////// RefCommon ////////////////////

namespace RefCommon {
	// Merge phased and non-phased VCF records
	// using phased positions as the reference.
	// For positions present in phased but missing in non-phased,
	// fill genotypes with "./.".
	// Positions not present in phased are ignored.
	std::vector<std::pair<std::size_t, std::vector<int>>>
	merge_records_core(const VCFGenoBase *phased_vcf,
						const VCFGenoBase *non_phased_vcf,
						const STRVEC& samples);
	
	std::vector<GenoRecord *> merge_records(const VCFGenoBase *phased_vcf,
											const VCFGenoBase *non_phased_vcf,
											const STRVEC& samples);
	std::vector<VCFFamilyRecord *>
	merge_family_records(const VCFGenoBase *phased_vcf,
							const VCFGenoBase *non_phased_vcf,
							const STRVEC& samples);
	
	// Merge phased and non-phased VCF records using phased positions
	// as reference.
	// Positions absent in the phased VCF are ignored.
	// Missing non-phased genotypes are filled with "./.".
    // Returns a list of (position, merged_genotypes).
   	std::vector<GenoRecord *> expand_records(const VCFGeno *vcf,
												const VCFGeno *phased_vcf);
};
#endif
