#ifndef __IMPUTEPROGONLY
#define __IMPUTEPROGONLY

#include <vector>
#include <string>
#include <utility>

#include "VCFImpFamilyRecord.h"

class TypeDeterminer;
class VCFHeteroHomoRecord;
class VCFImpFamilyRecord;
class VCFFamily;
class VCFHeteroHomoPP;
class VCFFillableRecord;
class Map;
class Option;


//////////////////// ImputeProgOnly ////////////////////

namespace ImputeProgOnly {
	void classify_record(std::size_t i, VCFFamilyRecord *record,
							const STRVEC& samples, const TypeDeterminer *td,
							std::vector<VCFHeteroHomoRecord *>& heho_records,
							std::vector<VCFImpFamilyRecord *>& other_records);
	std::pair<std::vector<VCFHeteroHomoRecord *>,
			  std::vector<VCFImpFamilyRecord *>>
	classify_records(VCFFamily *vcf, const STRVEC& samples, Option *option);
	VCFRecord *fill_NA(const VCFRecord *record1, const STRVEC& samples);
	VCFRecord *merge_record(const VCFRecord *record1,
										const VCFRecord *record2,
										const STRVEC& samples);
	VCFSmall *merge_parents_progenies(const VCFSmall *vcf_parents,
										const VCFSmall *vcf_progenies,
										const STRVEC& samples);
	VCFHeteroHomoPP *merge_vcf(
					std::array<std::vector<VCFFillableRecord *>, 4>& rss,
					const STRVEC& samples, const Map& gmap,
					const VCFSmall *vcf);
	VCFHeteroHomoPP *impute_prog_vcf_chr(const VCFSmall *parent_vcf,
										const VCFSmall *prog_vcf,
										const Map& gmap, const Option *Option);
}
#endif
