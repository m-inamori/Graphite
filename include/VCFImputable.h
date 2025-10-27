#ifndef __VCFIMPUTABLE
#define __VCFIMPUTABLE

#include "VCFGeno.h"
#include "Map.h"
#include "Haplotype.h"

class VCFSmall;


//////////////////// VCFImputable ////////////////////

class VCFImputable : public VCFGenoBase, public VCFMeasurable {
public:
	VCFImputable(const STRVEC& s, const Map& map_, const VCFSmall *vcf) :
									VCFGenoBase(s, vcf), VCFMeasurable(map_) { }
	virtual ~VCFImputable() { };
	
	virtual std::vector<Haplotype> collect_haplotypes_mat(
										std::size_t sample_index) const = 0;
	virtual std::vector<Haplotype> collect_haplotypes_pat(
										std::size_t sample_index) const = 0;
	
	std::vector<int> get_int_gts(std::size_t sample_index) const;
	Haplotype clip_haplotype(std::size_t sample_index, int side) const;
	bool is_block(const GenoRecord *record,
					std::vector<GenoRecord *>& rs) const {
		const double	length = cM(record->get_pos()) - cM(rs[0]->get_pos());
		if(length < 1.0)
			return true;
		else
			return rs.size() < 10 && length < 10.0;
	}
	
	HaplotypePair impute_cM_each_sample(HaplotypePair prev_hap,
										std::size_t sample_index, bool exec);
	void set_haplotype(HaplotypePair hap, std::size_t sample_index);
};
#endif
