#include <sstream>
#include <cassert>

#include "../include/VCFFillableRecord.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFFillableRecord ////////////////////

VCFFillableRecord *VCFFillableRecord::copy() const {
	return new VCFFillableRecord(pos, geno, index, type, comb, probs);
}

int VCFFillableRecord::gt_from_parent(int mat_from, int pat_from) const {
	if(mat_from == 0 || pat_from == 0)
		return Genotype::NA;
	
	const int	gt_from_mat = get_mat_allele(mat_from-1);
	const int	gt_from_pat = get_pat_allele(pat_from-1);
	return gt_from_mat | (gt_from_pat << 1) | 4;
}

int VCFFillableRecord::gt_from_mat(int mat_from, size_t i) const {
	const int	prog_gt = unphased(i);
	// this allele must come from mat
	const int	mat_allele = get_allele(0, mat_from-1);
	const int	pat_gt = this->pat_gt();
	if(prog_gt == Genotype::UN_00) {
		if(mat_allele == 0)
			return Genotype::PH_00;
		else if(pat_gt < Genotype::PH_11)
			return Genotype::PH_10;
		else
			return Genotype::PH_11;
	}
	else if(prog_gt == Genotype::UN_11) {
		if(mat_allele == 1)
			return Genotype::PH_11;
		else if(pat_gt == Genotype::PH_00)
			return Genotype::PH_00;
		else
			return Genotype::PH_01;
	}
	else if(prog_gt == Genotype::UN_01) {
		if(mat_allele == 0)
			return Genotype::PH_01;
		else
			return Genotype::PH_10;
	}
	else {
		return Genotype::NA;
	}
}

int VCFFillableRecord::gt_from_pat(int pat_from, size_t i) const {
	const int	prog_gt = unphased(i);
	const int	pat_allele = get_allele(1, pat_from-1);
	const int	mat_gt = this->mat_gt();
	if(prog_gt == Genotype::UN_00) {
		if(pat_allele == 0)
			return Genotype::PH_00;
		else if(mat_gt == Genotype::PH_11)
			return Genotype::PH_11;
		else
			return Genotype::PH_01;
	}
	else if(prog_gt == Genotype::UN_11) {
		if(pat_allele == 1)
			return Genotype::PH_11;
		else if(mat_gt == Genotype::PH_00)
			return Genotype::PH_10;
		else
			return Genotype::PH_00;
	}
	else if(prog_gt == Genotype::UN_01) {
		if(pat_allele == 0)
			return Genotype::PH_10;
		else
			return Genotype::PH_01;
	}
	else {
		return Genotype::NA;
	}
}

int VCFFillableRecord::mat_from(size_t i) const {
	if(this->is_NA(i))
		return 0;
	else if(!this->is_mat_hetero())
		return 0;
	else {
		const int	a = this->get_mat_allele(0);
		// TODO: if the genotype is non-phased, something wrong.
		const int	gt = this->geno[i];
		const int	b = (gt == Genotype::UN_11 ||
						 gt == Genotype::PH_10 ||
						 gt == Genotype::PH_11) ? 1: 0;
		
		if(a == b)
			return 1;
		else
			return 2;
	}
}

int VCFFillableRecord::pat_from(size_t i) const {
	if(this->is_NA(i))
		return 0;
	else if(!this->is_pat_hetero())
		return 0;
	else {
		const int	a = this->get_pat_allele(0);
		// TODO: if the genotype is non-phased, something wrong.
		const int	gt = this->geno[i];
		const int	b = (gt == Genotype::UN_00 ||
						 gt == Genotype::PH_00 ||
						 gt == Genotype::PH_10) ? 0 : 1;
		
		if(a == b)
			return 1;
		else
			return 2;
	}
}

void VCFFillableRecord::inverse_parents_gts(bool inv_mat, bool inv_pat) {
	// both must be hetero
	if(inv_mat)
		this->geno[0] = Genotype::inverse(this->geno[0]);
	if(inv_pat)
		this->geno[1] = Genotype::inverse(this->geno[1]);
}

bool VCFFillableRecord::is_near_prog_gts(const vector<int>& gts) const {
	int	num = 0;
	int	dist = 0;
	for(size_t i = 0; i < gts.size(); ++i) {
		// TODO: modify
		if(geno[i+2] != Genotype::UN_01)
			num += 1;
		if(geno[i+2] == Genotype::UN_01) {
			if(!(gts[i] == Genotype::PH_01 || gts[i] == Genotype::PH_10))
				dist += 1;
		}
		else if(geno[i+2] == Genotype::UN_00) {
			if(gts[i] != Genotype::PH_00)
				dist += 1;
		}
		else if(geno[i+2] == Genotype::UN_11) {
			if(gts[i] != Genotype::PH_11)
				dist += 1;
		}
		else {
			dist += 1;
		}
		/*
		if(Genotype::is_NA(geno[i+2])
			dist += 1;
		else if(Genotype::unphased(gts[i]) != geno[i+2])
			dist += 1;
		*/
	}
	return dist < max(1, num / 2);
}

int VCFFillableRecord::distance(const vector<int>& gts1,
								const vector<int>& gts2) {
	int	dist = 0;
	for(size_t i = 0; i < gts1.size(); ++i) {
		if(!Genotype::is_same_gts(gts1[i], gts2[i]))
			dist += 1;
	}
	return dist;
}

void VCFFillableRecord::modify_gts(const vector<int>& new_prog_gts) {
	vector<int>	orig_gts(geno.begin() + 2, geno.end());
	int	min_dist = distance(new_prog_gts, orig_gts);
	vector<int>	min_gts = new_prog_gts;
	int	min_i = 0;
	for(int i = 1; i < 4; ++i) {
		const bool	inv_mat = (i & 2) == 2;
		const bool	inv_pat = (i & 1) == 1;
		const vector<int>	inv_prog_gts = inverse_prog_gts(new_prog_gts,
															inv_mat, inv_pat);
		const int	dist = distance(inv_prog_gts, orig_gts);
		if(dist < min_dist) {
			min_dist = dist;
			min_gts = inv_prog_gts;
			min_i = i;
		}
	}
	
	const bool	inv_mat = (min_i & 2) == 2;
	const bool	inv_pat = (min_i & 1) == 1;
	this->inverse_parents_gts(inv_mat, inv_pat);
	std::copy(min_gts.begin(), min_gts.end(), this->geno.begin() + 2);
}

int VCFFillableRecord::inverse_prog_gt(int gt,
										bool inv_mat, bool inv_pat) {
	if(Genotype::is_NA(gt))
		return Genotype::NA;
	
	if(gt == 0)
		gt = 4;
	else if(gt == 1)
		gt = 6;
	else if(gt == 2)
		gt = 7;
	
	const int	a1 = Genotype::get_allele(gt, 0);
	const int	a2 = Genotype::get_allele(gt, 1);
	return Genotype::from_alleles(Genotype::inverse_allele(a1, inv_mat),
								  Genotype::inverse_allele(a2, inv_pat));
}

vector<int> VCFFillableRecord::inverse_prog_gts(const vector<int>& prog_gts,
											bool inv_mat, bool inv_pat) const {
	vector<int>	inv_prog_gts;
	for(auto p = prog_gts.begin(); p != prog_gts.end(); ++p) {
		inv_prog_gts.push_back(inverse_prog_gt(*p, inv_mat, inv_pat));
	}
	return inv_prog_gts;
}

void VCFFillableRecord::modify_parents_type() {
	const int	mat_gt = this->mat_gt();
	const int	pat_gt = this->pat_gt();
	if(this->comb != ParentComb::P00x11 &&
			((mat_gt == Genotype::PH_00 && pat_gt == Genotype::PH_11) ||
			 (mat_gt == Genotype::PH_11 && pat_gt == Genotype::PH_00)))
		comb = ParentComb::P00x11;
}

int VCFFillableRecord::from_which_chrom(size_t i, bool is_mat) const {
	const int j = is_mat ? 0 : 1;
	const int	parent_allele = this->get_allele(j, 0);
	return parent_allele == this->get_allele(i, j) ? 1 : 2;
}

int VCFFillableRecord::hash(int d) const {
	int	hash_ = 0;
	int	n = this->get_pos();
	while(n > 0) {
		const int	r = n % 10;
		n /= 10;
		hash_ += r;
	}
	return hash_ % d;
}

int VCFFillableRecord::decide_by_majority(const vector<int>& GTs) const {
	// ex) [0/0, 0/1, 0/0, 1/1] -> { 0/0: 2, 0/1: 1, 1/1: 1 }
	map<int, int>	counter;
	for(auto p = GTs.begin(); p != GTs.end(); ++p)
		counter[*p] += 1;
	
	// ex) { 0/0: 2, 0/1: 1, 1/1: 1 } -> { 2: [0/0], 1: [0/1, 1/1] }
	map<int, vector<int>>	dic_num;
	for(auto p = counter.begin(); p != counter.end(); ++p)
		dic_num[p->second].push_back(p->first);
	
	pair<int, vector<int>>	max_pair;
	for(auto p = dic_num.begin(); p != dic_num.end(); ++p) {
		if(p->first > max_pair.first)
			max_pair = *p;
	}
	
	vector<int>	max_GTs = max_pair.second;
	std::sort(max_GTs.begin(), max_GTs.end());	// to match the Python version
	if(max_GTs.size() == 1)
		return max_GTs.front();
	else	// Use a random number if there is more than one candidates
		return max_GTs[this->hash((int)max_GTs.size())];
}

void VCFFillableRecord::swap_parents(size_t i, int GT) {
	if(GT == this->geno[i])
		return;
	
	if(GT == Genotype::PH_00 || GT == Genotype::PH_11) {
		const bool	is_mat_00 = (i == 0) ^ (GT == Genotype::PH_11);
		this->geno[0] = is_mat_00 ? Genotype::PH_00 : Genotype::PH_11;
		this->geno[1] = is_mat_00 ? Genotype::PH_11 : Genotype::PH_00;
		
		// swap genotyeps of progenies
		const int	prog_GT = is_mat_00 ? Genotype::PH_01 : Genotype::PH_10;
		for(size_t j = 2; j < this->geno.size(); ++j)
			this->geno[j] = prog_GT;
	}
}

// Determine a single Genotype
// when the Genotype of the same sample differs from family to family
int VCFFillableRecord::decide_duplicated_Genotype(
								const vector<VCFFillableRecord *>& records,
								const vector<pair<size_t, size_t>>& positions) {
	vector<int>	GTs;
	for(auto p = positions.begin(); p != positions.end(); ++p)
		GTs.push_back(records[p->first]->geno[p->second]);
	
	// if all genotype are ./., not change
	if(std::all_of(GTs.begin(), GTs.end(),
			static_cast<bool(*)(int)>(&Genotype::is_NA)))
		return Genotype::NA;
	
	// if the sample is a progeny of a family, use the genotype
	for(size_t i = 0; i < GTs.size(); ++i) {
		const auto	j = positions[i].second;	// jth(0-based) sample of family
		const int	GT = GTs[i];
		if(j >= 2 && !Genotype::is_NA(GT))		// j >= 2 means progeny
			return GT;
	}
	
	vector<int>	GTs_less_NA;
	for(auto p = GTs.begin(); p != GTs.end(); ++p) {
		if(!Genotype::is_NA(*p))
			GTs_less_NA.push_back(*p);
	}
	
	// If all Genotypes except ./. are the same Genotype, use that Genotype
	if(Common::is_all_same(GTs_less_NA))
		return GTs_less_NA.front();
	
	vector<int>	GTs_less_00x11;
	for(size_t i = 0; i < GTs.size(); ++i) {
		const auto	*record = records[positions[i].first];
		const int	GT = GTs[i];
		if(!Genotype::is_NA(GT) && record->comb != ParentComb::P00x11)
			GTs_less_00x11.push_back(GT);
	}
	
	if(GTs_less_00x11.empty())	// only 0/0 x 1/1
		return records.front()->decide_by_majority(GTs_less_NA);
	else if(Common::is_all_same(GTs_less_00x11))
		return GTs_less_00x11.front();
	else
		return records.front()->decide_by_majority(GTs_less_00x11);
}

vector<VCFFillableRecord *> VCFFillableRecord::convert(
									const vector<VCFImpFamilyRecord *>& records,
									const VCFHeteroHomo *vcf) {
	const VCFSmall	*ref_vcf = vcf->get_ref_vcf();
	const auto	cols = ref_vcf->extract_columns(vcf->get_samples());
	vector<VCFFillableRecord *>	new_records;
	for(size_t i = 0; i < records.size(); ++i) {
		const auto	*record = records[i];
		const auto	type = record->get_fill_type();
		const auto	*ref_record = ref_vcf->get_record(i);
		const auto	probs = ref_record->parse_PL(record->get_geno(), cols);
		auto	*r = new VCFFillableRecord(record->get_pos(), record->get_geno(),
											record->get_index(), type,
											record->get_comb(), probs);
		new_records.push_back(r);
	}
	return new_records;
}

GenoRecord *VCFFillableRecord::merge(const vector<VCFFillableRecord *>& records,
														const STRVEC& samples) {
	const ll	pos = records[0]->get_pos();
	vector<int>	new_geno;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const vector<int>&	geno = (*p)->get_geno();
		new_geno.insert(new_geno.end(), geno.begin(), geno.end());
	}
	return new GenoRecord(pos, new_geno);
}

void VCFFillableRecord::integrate_each_sample(
								const vector<VCFFillableRecord *>& records,
								const vector<pair<size_t, size_t>>& positions) {
	const int	GT = decide_duplicated_Genotype(records, positions);
	for(auto p = positions.begin(); p != positions.end(); ++p) {
		VCFFillableRecord	*record = records[p->first];
		if(record->is_00x11() && p->second <= 1)
			record->swap_parents(p->second, GT);
	}
}

bool VCFFillableRecord::is_all_same_GT(
						const vector<VCFFillableRecord *>& records,
						const vector<pair<size_t, size_t>>& pos_samples) {
	if(pos_samples.size() == 1)
		return true;
	
	const auto	pos_sample = pos_samples.front();
	const int	GT0 = records[pos_sample.first]->geno[pos_sample.second];
	for(auto p = pos_samples.begin() + 1; p != pos_samples.end(); ++p) {
		const int	GT = records[p->first]->geno[p->second];
		if(GT != GT0)
			return false;
	}
	return true;
}

GenoRecord *VCFFillableRecord::integrate(
					const vector<VCFFillableRecord *>& records,
					const vector<string>& samples,
					const vector<vector<pair<size_t, size_t>>>& pos_samples) {
	// Exchangeable if one parent Genotype is 0|0 and the other is 1|1
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		if(!is_all_same_GT(records, *p)) {
			integrate_each_sample(records, *p);
		}
	}
	
	const ll	pos = records[0]->get_pos();
	vector<int>	geno;
	for(auto p = pos_samples.begin(); p != pos_samples.end(); ++p) {
		const pair<int, int>&	pos = p->front();
		geno.push_back(records[pos.first]->geno[pos.second]);
	}
	return new GenoRecord(pos, geno);
}

// after imputed
int VCFFillableRecord::from_which_chrom(const VCFFillableRecord *record,
													size_t i, bool is_mat) {
	if(record == NULL)
		return 0;
	
	return record->from_which_chrom(i, is_mat);
}
