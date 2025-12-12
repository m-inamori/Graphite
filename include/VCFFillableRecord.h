#ifndef __VCFFILLABLERECORD
#define __VCFFILLABLERECORD

#include <vector>
#include <string>

#include "VCFImpFamilyRecord.h"
#include "ClassifyRecord.h"

class VCFHeteroHomo;


//////////////////// VCFFillableRecord ////////////////////

class VCFFillableRecord : public VCFFamilyRecord {
protected:
	int			index;
	FillType	type;
	ParentComb	comb;
	std::vector<VCFRecord::Probs>	probs;
	
public:
	VCFFillableRecord(ll pos, const std::vector<int>& geno,
						int i, FillType t, ParentComb c,
						const std::vector<VCFRecord::Probs>& ps) :
								VCFFamilyRecord(pos, geno), index(i),
								type(t), comb(c), probs(ps) { }
	~VCFFillableRecord() { }
	
	int get_index() const { return index; }
	FillType get_type() const { return type; }
	bool is_unable() const { return type == FillType::UNABLE; }
	bool is_fillable_type() const {
		return type == FillType::UNABLE || type == FillType::IMPUTABLE;
	}
	bool is_filled_type() const {
		return type == FillType::FILLED || type == FillType::MAT
										|| type == FillType::PAT;
	}
	bool is_mat_type() const { return type == FillType::MAT; }
	bool is_pat_type() const { return type == FillType::PAT; }
	
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
	
	VCFFillableRecord *copy() const;
	int gt_from_parent(int mat_from, int pat_from) const;
	int gt_from_mat(int mat_from, std::size_t i) const;
	int gt_from_pat(int pat_from, std::size_t i) const;
	int mat_from(std::size_t i) const;
	int pat_from(std::size_t i) const;
	
	void modify_parents_type();
	int from_which_chrom(std::size_t i, bool is_mat) const;
	void modify_gts(const std::vector<int>& new_prog_gts);
	
private:
	std::vector<int> inverse_prog_gts(const std::vector<int>& prog_gts,
											bool inv_mat, bool inv_pat) const;
	void inverse_parents_gts(bool inv_mat, bool inv_pat);
	bool is_near_prog_gts(const std::vector<int>& gts) const;
	int hash(int d) const;
	int decide_by_majority(const std::vector<int>& GTs) const;
	void swap_parents(std::size_t i, int GT);
	
	static int distance(const std::vector<int>& gts1,
								const std::vector<int>& gts2);
	static int inverse_prog_gt(int gt, bool inv_mat, bool inv_pat);
	
public:
	static int decide_duplicated_Genotype(
			const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::pair<std::size_t, std::size_t>>& positions);
	GenoRecord *merge(const std::vector<VCFFillableRecord *>& records,
														const STRVEC& samples);
	static std::vector<VCFFillableRecord *> convert(
							const std::vector<VCFImpFamilyRecord *>& records,
							const VCFHeteroHomo *vcf);
	static void integrate_each_sample(
			const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::pair<std::size_t, std::size_t>>& positions);
	static bool is_all_same_GT(const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::pair<std::size_t, std::size_t>>& pos_samples);
	static GenoRecord *integrate(const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::string>& samples,
			const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& pos_samples);
	static int from_which_chrom(const VCFFillableRecord *record,
										std::size_t i, bool is_mat);
	static int from_which_chrom_mat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, true);
	}
	static int from_which_chrom_pat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, false);
	}
};
#endif
