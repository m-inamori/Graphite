#ifndef __VCFSELFFILLABLERECORD
#define __VCFSELFFILLABLERECORD

#include <vector>
#include <string>

#include "VCF.h"
#include "VCFImpSelfRecord.h"
#include "ClassifyRecord.h"


//////////////////// VCFSelfFillableRecord ////////////////////

class VCFSelfFillableRecord : public GenoRecord {
protected:
	int			index;
	SelfFillType	type;
	ParentComb	comb;
	std::vector<VCFRecord::Probs>	probs;
	
public:
	VCFSelfFillableRecord(ll pos, const std::vector<int>& geno,
							int i, SelfFillType t, ParentComb c,
							const std::vector<VCFRecord::Probs>& ps) :
				GenoRecord(pos, geno), index(i), type(t), comb(c), probs(ps) { }
	~VCFSelfFillableRecord() { }
	
	int get_index() const { return index; }
	SelfFillType get_type() const { return type; }
	bool is_unable() const { return type == SelfFillType::UNABLE; }
	bool is_fillable_type() const {
		return type == SelfFillType::UNABLE || type == SelfFillType::IMPUTABLE;
	}
	bool is_filled_type() const {
		return type == SelfFillType::FILLED || type == SelfFillType::P01;
	}
	
	ParentComb get_comb() const { return comb; }
	bool is_00x11() const { return comb == ParentComb::P00x11; }
	void set_00x11() { comb = ParentComb::P00x11; }
	
	double get_prob(std::size_t i, int gt) const {
		switch(gt) {
			case 0:  return std::get<0>(probs[i]);
			case 1:  return std::get<1>(probs[i]);
			default: return std::get<2>(probs[i]);
		}
	}
	
	int from_which_chrom(std::size_t i, bool is_mat) const;
	
public:
	static VCFSelfFillableRecord *convert(const VCFImpSelfRecord *record,
									const std::vector<VCFRecord::Probs>& probs);
	static int from_which_chrom(const VCFSelfFillableRecord *record,
										std::size_t i, bool is_mat);
	static int from_which_chrom_mat(const VCFSelfFillableRecord *record,
															std::size_t i) {
		return from_which_chrom(record, i, true);
	}
	static int from_which_chrom_pat(const VCFSelfFillableRecord *record,
															std::size_t i) {
		return from_which_chrom(record, i, false);
	}
};
#endif
