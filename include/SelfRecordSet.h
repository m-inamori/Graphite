// SelfRecordSet.h
#ifndef __SELFRECORDSET
#define __SELFRECORDSET

#include <vector>
#include <string>

#include "VCFSelfFillableRecord.h"
#include "TypeDeterminer.h"
#include "Genotype.h"


//////////////////// SelfRecordSet ////////////////////

class SelfRecordSet {
public:
	using Pair = std::pair<int,int>;
	
public:
	VCFSelfFillableRecord	*record;
	VCFSelfFillableRecord	*prev_record;
	VCFSelfFillableRecord	*next_record;
	
public:
	SelfRecordSet(VCFSelfFillableRecord *r,
			  VCFSelfFillableRecord *p, VCFSelfFillableRecord *n) :
				  				record(r), prev_record(p), next_record(n) { }
	virtual ~SelfRecordSet() { }
	
	std::string gt_each(int i, const VCFSelfFillableRecord *record) const {
		return record == NULL ? "" : record->get_gt(i);
	}
	std::string gt(int i) const { return gt_each(i, record); }
	std::string prev_mat_gt(int i) const { return gt_each(i, prev_record); }
	std::string next_mat_gt(int i) const { return gt_each(i, next_record); }
	std::string prev_pat_gt(int i) const { return gt_each(i, prev_record); }
	std::string next_pat_gt(int i) const { return gt_each(i, next_record); }
	
	int prev_mat_from(std::size_t i) const;
	int next_mat_from(std::size_t i) const;
	int prev_pat_from(std::size_t i) const;
	int next_pat_from(std::size_t i) const;
	
	int near_mat_from(size_t i) const;
	int near_pat_from(size_t i) const;
	
	bool is_prev_near() const;
	
	// phasingされている前提
	int from_which_chrom(const std::string& gt,
							const VCFRecord *record, bool mat) const {
		if(record == NULL)
			return 0;
		
		const std::size_t	i = mat ? 0 : 1;
		const std::string&	parent_gt = record->get_gt(i);
		return parent_gt.c_str()[0] == gt[i*2] ? 1 : 2;
	}
	
	int from_which_chrom_prev_mat(const std::string& gt) const {
		return from_which_chrom(gt, prev_record, true);
	}
	int from_which_chrom_next_mat(const std::string& gt) const {
		return from_which_chrom(gt, next_record, true);
	}
	int from_which_chrom_prev_pat(const std::string& gt) const {
		return from_which_chrom(gt, prev_record, false);
	}
	int from_which_chrom_next_pat(const std::string& gt) const {
		return from_which_chrom(gt, next_record, false);
	}
	
	std::vector<double> likelihoods_from_which_chrom(int prev_chrom,
												int next_chrom) const;
	std::vector<double> likelihoods_from_which_chrom(std::size_t i,
												bool is_mat) const;
	double compute_phasing_likelihood_each(int phasing, std::size_t i) const;
	double likelihood_each(int phasing, std::size_t i,
							const std::vector<double>& probs_mat,
							const std::vector<double>& probs_pat) const;
	double compute_phasing_likelihood(int phasing) const;
	int select_phasing(const std::vector<int>& candidates) const;
	int determine_phasing_core(
					const std::vector<std::pair<double, int>>& lls) const;
	std::vector<int> possible_phasings() const;
	void determine_parent_phasing() const;
};
#endif
