// RecordSet.h
#ifndef __RECORDSET
#define __RECORDSET

#include <vector>
#include <string>

#include "VCFFillableRecord.h"
#include "TypeDeterminer.h"
#include "Genotype.h"


//////////////////// RecordSet ////////////////////

class RecordSet {
public:
	using Pair = std::pair<int,int>;
	
public:
	VCFFillableRecord	*record;
	VCFFillableRecord	*prev_mat_record;
	VCFFillableRecord	*next_mat_record;
	VCFFillableRecord	*prev_pat_record;
	VCFFillableRecord	*next_pat_record;
	
public:
	RecordSet(VCFFillableRecord *r,
			  VCFFillableRecord *pm, VCFFillableRecord *nm,
			  VCFFillableRecord *pp, VCFFillableRecord *np) : record(r),
					  			prev_mat_record(pm), next_mat_record(nm),
					  			prev_pat_record(pp), next_pat_record(np) { }
	virtual ~RecordSet() { }
	
	int gt_each(std::size_t i, const VCFFillableRecord *record) const {
		return record == NULL ? Genotype::NA : record->get_geno()[i];
	}
	int gt(std::size_t i) const { return gt_each(i, record); }
	int prev_mat_gt(std::size_t i) const { return gt_each(i, prev_mat_record); }
	int next_mat_gt(std::size_t i) const { return gt_each(i, next_mat_record); }
	int prev_pat_gt(std::size_t i) const { return gt_each(i, prev_pat_record); }
	int next_pat_gt(std::size_t i) const { return gt_each(i, next_pat_record); }
	
	int prev_mat_from(std::size_t i) const;
	int next_mat_from(std::size_t i) const;
	int prev_pat_from(std::size_t i) const;
	int next_pat_from(std::size_t i) const;
	
	int near_mat_from(size_t i) const;
	int near_pat_from(size_t i) const;
	
	bool is_mat_prev_near() const;
	bool is_pat_prev_near() const;
	bool is_prev_nearer(bool is_mat) const;
	
	// phasingされている前提
	int from_which_chrom(int gt, const VCFFillableRecord *record,
															bool mat) const {
		if(record == NULL)
			return 0;
		
		const std::size_t	i = mat ? 0 : 1;
		const int	parent_gt = record->get_geno()[i];
		return (parent_gt & 1) == ((gt >> i) & 1) ? 1 : 2;
	}
	
	int from_which_chrom_prev_mat(int gt) const {
		return from_which_chrom(gt, prev_mat_record, true);
	}
	int from_which_chrom_next_mat(int gt) const {
		return from_which_chrom(gt, next_mat_record, true);
	}
	int from_which_chrom_prev_pat(int gt) const {
		return from_which_chrom(gt, prev_pat_record, false);
	}
	int from_which_chrom_next_pat(int gt) const {
		return from_which_chrom(gt, next_pat_record, false);
	}
	
	Pair select_nearest_froms(const std::vector<Pair>& pairs,
												std::size_t i) const;
	Pair select_pair(const std::vector<Pair>& pairs,
							std::size_t i, bool selected=false) const;
	std::vector<double> likelihoods_from_which_chrom(int prev_chrom,
												int next_chrom) const;
	std::vector<double> likelihoods_from_which_chrom(std::size_t i,
												bool is_mat) const;
	double compute_phasing_likelihood_each(int mat_phasing,
										int pat_phasing, std::size_t i) const;
	double likelihood_each(const std::vector<double>& probs_mat,
							const std::vector<double>& probs_pat,
							int mat_phasing, int pat_phasing,
							std::size_t i) const;
	virtual double compute_phasing_likelihood(int mat_phasing,
													int pat_phasing) const;
	std::pair<int, int> select_phasing(
					const std::vector<std::pair<int, int>>& candidates) const;
	std::pair<int, int> determine_phasing_core(
					const std::vector<std::tuple<double, int, int>>& lls) const;
	void determine_phasing() const;
	virtual std::vector<std::pair<int, int>> possible_phasings() const;
	int select_from(int from1, int from2, const GenoRecord *record1,
											const GenoRecord *record2) const;
	int modify_gt(size_t i) const;
	void impute(bool necessary_parents_phasing) const;
	void impute_core() const;
	
	int select_mat(const std::vector<Pair>& pairs) const;
	void impute_NA_mat_each(std::size_t i) const;
	int select_pat(const std::vector<Pair>& pairs) const;
	void impute_NA_pat_each(std::size_t i) const;
	int determine_mat_from(std::size_t i) const;
	int determine_pat_from(std::size_t i) const;
};


//////////////////// RecordSetSmall ////////////////////

class RecordSetSmall : public RecordSet {
public:
	RecordSetSmall(VCFFillableRecord *r,
				   VCFFillableRecord *pm, VCFFillableRecord *nm,
				   VCFFillableRecord *pp, VCFFillableRecord *np) :
								   	RecordSet(r, pm, nm, pp, np) { }
	~RecordSetSmall() { }
	
	std::vector<std::pair<int, int>> possible_phasings() const override;
	double compute_phasing_likelihood(int mat_phasing,
									  int pat_phasing) const override;
};
#endif
