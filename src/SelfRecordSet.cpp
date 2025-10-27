#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../include/SelfRecordSet.h"
#include "../include/log.h"
#include "../include/common.h"

using namespace std;


//////////////////// SelfRecordSet ////////////////////

int SelfRecordSet::prev_mat_from(std::size_t i) const {
	if(prev_record == NULL)
		return 0;
	return prev_record->from_which_chrom(i, true);
}
int SelfRecordSet::next_mat_from(std::size_t i) const {
	if(next_record == NULL)
		return 0;
	return next_record->from_which_chrom(i, true);
}
int SelfRecordSet::prev_pat_from(std::size_t i) const {
	if(prev_record == NULL)
		return 0;
	return prev_record->from_which_chrom(i, false);
}
int SelfRecordSet::next_pat_from(std::size_t i) const {
	if(next_record == NULL)
		return 0;
	return next_record->from_which_chrom(i, false);
}

// Likelihood based on which of the parent Haplotypes
// came from before and after the record
vector<double> SelfRecordSet::likelihoods_from_which_chrom(
										int prev_from, int next_from) const {
	static const double ps[] = {
		0.5, 0.9, 0.1, 0.9, 0.99, 0.5, 0.1, 0.5, 0.01
	};
	const double	prob1 = ps[prev_from+next_from*3];
	vector<double>	probs = { prob1, 1.0 - prob1 };
	return probs;
}

vector<double> SelfRecordSet::likelihoods_from_which_chrom(
												size_t i, bool is_mat) const {
	if(is_mat)
		return likelihoods_from_which_chrom(prev_mat_from(i), next_mat_from(i));
	else
		return likelihoods_from_which_chrom(prev_pat_from(i), next_pat_from(i));
}

double SelfRecordSet::likelihood_each(int phasing, size_t i,
										const vector<double>& probs_mat,
										const vector<double>& probs_pat) const {
	double	likelihood = 0.0;
	for(int l = 0; l < 4; ++l) {
		const int	j = l >> 1;
		const int	k = l & 1;
		const int	gt = ((phasing >> j) & 1) + ((phasing >> k) & 1);
		likelihood += probs_mat[j] * probs_pat[k] * record->get_prob(i, gt);
	}
	return Log::modified_log(likelihood);
}

double SelfRecordSet::compute_phasing_likelihood_each(int phasing,
															size_t i) const {
	const auto	probs_mat = this->likelihoods_from_which_chrom(i, true);
	const auto	probs_pat = this->likelihoods_from_which_chrom(i, false);
	return this->likelihood_each(phasing, i, probs_mat, probs_pat);
}

double SelfRecordSet::compute_phasing_likelihood(int phasing) const {
	if(this->record == NULL)
		return log(0.0001);
	
	double	ll = 0.0;
	for(int i = 1; i < (int)record->num_samples(); ++i) {
		ll += compute_phasing_likelihood_each(phasing, i);
	}
	return ll;
}

int SelfRecordSet::select_phasing(const vector<int>& candidates) const {
	if(candidates.size() == 1 || this->record == NULL)
		return candidates[0];
	
	// Select the genotype closest to the original
	const int	parent_gt = this->record->unphased(0);
	vector<pair<int, int>>	v;
	for(auto p = candidates.begin(); p != candidates.end(); ++p) {
		const int	phasing = *p;
		const int	parent_int_gt1 = (phasing >> 1) + (phasing & 1);
		const int	score = abs(parent_int_gt1 - parent_gt);
		v.push_back(make_pair(score, phasing));
	}
	
	auto	p = min_element(v.begin(), v.end());
	return p->second;
}

int SelfRecordSet::determine_phasing_core(
							const vector<pair<double, int>>& lls) const {
	// If there are almost identical values, collect them and select again.
	double	max_ll = lls.back().first;
	vector<int>	candidates(1, lls.back().second);
	for(size_t j = 1; j < lls.size(); ++j) {
		const size_t	i = lls.size() - 1 - j;
		if(lls[i].first > max_ll - 1e-9)
			candidates.push_back(lls[i].second);	// phasing
		else
			break;
	}
	return select_phasing(candidates);
}

vector<int> SelfRecordSet::possible_phasings() const {
	if(record == NULL)
		return vector<int>();
	
	switch(this->record->get_comb()) {
		case ParentComb::P00x00: return { 0 };
		case ParentComb::P11x11: return { 3 };
		default:                 return { 0, 1, 2, 3 };
	}
}

void SelfRecordSet::determine_parent_phasing() const {
	if(this->record == NULL)
		return;
	
	const vector<int>	phasings = this->possible_phasings();
	vector<pair<double, int>>	lls;
	for(auto p = phasings.begin(); p != phasings.end(); ++p) {
		const double	ll = compute_phasing_likelihood(*p);
		lls.push_back(make_pair(ll, *p));
	}
	
	// Measures for when there is almost the same ll as the largest ll
	std::sort(lls.begin(), lls.end());
	const int	phasing = this->determine_phasing_core(lls);
	this->record->set_geno(0, phasing | 4);
}
