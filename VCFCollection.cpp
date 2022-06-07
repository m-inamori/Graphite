#include <random>
#include <cassert>
#include "VCFCollection.h"
#include "option.h"
#include "Map.h"
#include "BiasProbability.h"
#include "common.h"

using namespace std;


//////////////////// VCFCollection ////////////////////

VCFCollection::VCFCollection(const vector<VCFImputable *>& vcfs_) :
														vcfs(vcfs_) {
	for(int i = 0; i < (int)vcfs.size(); ++i)
		vcfs[i]->set_group_id(i);
}

VCFCollection::~VCFCollection() {
	// VCFHeteroHomo::divide_into_chromosomesで作ったMapはvcfsで共用
	// 一つ消せばよい
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		delete *p;
}

size_t VCFCollection::max_index() const {
	size_t	index = vcfs.front()->max_index();
	for(auto p = vcfs.begin() + 1; p != vcfs.end(); ++p) {
		index = std::max((*p)->max_index(), index);
	}
	return index;
}

void VCFCollection::impute(int min_positions, double MIN_CROSSOVER, int T) {
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		if((int)(*p)->size() >= min_positions)
			(*p)->impute(MIN_CROSSOVER, T);
	}
}

// topだけはfalse
vector<vector<bool>> VCFCollection::all_boolean_vectors(size_t n, bool top) {
	// D&C
	if(n == 1U) {
		if(top) {
			vector<vector<bool>>	bs = { vector<bool>(1, false) };
			return bs;
		}
		else {
			vector<vector<bool>>	bs = { vector<bool>(1, true),
											vector<bool>(1, false) };
			return bs;
		}
	}
	else if(n % 2 == 1U) {
		const auto	part_vectors = all_boolean_vectors(n - 1, top);
		vector<vector<bool>>	vectors;
		for(auto p = part_vectors.begin(); p != part_vectors.end(); ++p) {
			vector<bool>	v = *p;
			v.push_back(true);
			vectors.push_back(v);
			v.pop_back();
			v.push_back(false);
			vectors.push_back(v);
		}
		return vectors;
	}
	else {
		const auto	part_vectors1 = all_boolean_vectors(n / 2, top);
		const auto	part_vectors2 = all_boolean_vectors(n / 2, false);
		vector<vector<bool>>	vectors;
		for(auto p = part_vectors1.begin(); p != part_vectors1.end(); ++p) {
			vector<bool>	v1 = *p;
			for(auto q = part_vectors2.begin(); q != part_vectors2.end(); ++q) {
				vector<bool>	v2 = *q;
				v1.insert(v1.end(), v2.begin(), v2.end());
				vectors.push_back(v1);
				v1.erase(v1.begin() + n/2, v1.end());
			}
		}
		return vectors;
	}
}

vector<vector<bool>> VCFCollection::make_random_boolean_vectors(
												size_t num_vec, size_t n) {
	std::mt19937	mt(123456789);
	std::uniform_int_distribution<int>	dist(0, 1);
	
	vector<vector<bool>>	vectors(num_vec, vector<bool>(n));
	for(size_t i = 0; i < num_vec; ++i) {
		vectors[i][0] = false;
		for(size_t j = 1; j < n; ++j)
			vectors[i][j] = dist(mt) == 0;
	}
	return vectors;
}

vector<vector<bool>> VCFCollection::make_boolean_vectors(size_t n) {
	const size_t	bs_limit = 10;
	if(vcfs.size() <= bs_limit) {
		return all_boolean_vectors(vcfs.size(), true);
	}
	else {
		return make_random_boolean_vectors(1U << bs_limit, vcfs.size());
	}
}

// 元のrecordはどうVCFに分配されているのか
// [(index of vcfs, index in vcf)]
vector<pair<size_t,size_t>> VCFCollection::make_which_vcf_table() const {
	const size_t	num_markers = max_index() + 1;
	vector<pair<size_t,size_t>>	which_vcfs(num_markers);
	for(size_t i = 0U; i < vcfs.size(); ++i) {
		const auto&	vcf = vcfs[i];
		for(size_t j = 0U; j < vcf->size(); ++j) {
			const auto	*record = vcf->get_record(j);
			which_vcfs[record->get_index()] = pair<size_t,size_t>(i, j);
		}
	}
	return which_vcfs;
}

vector<VCFCollection::Joint> VCFCollection::extract_joints() const {
	vector<Joint>	joints;
	const auto	which_vcfs = make_which_vcf_table();
	for(auto p = which_vcfs.begin(); p != which_vcfs.end() - 1; ++p) {
		const auto&	pair1 = *p;
		const auto&	pair2 = *(p + 1);
		if(pair1.first != pair2.first)	// different VCF indices
			joints.push_back(Joint(pair1, pair2));
	}
	return joints;
}

// ヘテロ親のハプロタイプを逆にしないときまたは逆にしたときの繋がり具合
// 各サンプルで繋がっていないと+1
int VCFCollection::score_joint(const Joint& joint, bool b) const {
	// b: trueが逆にしないとき
	const auto&	pair1 = joint.first;
	const auto&	pair2 = joint.second;
	auto	*record1 = vcfs[pair1.first]->get_record(pair1.second);
	auto	*record2 = vcfs[pair2.first]->get_record(pair2.second);
	return record1->difference(record2, !b);
}

// VCF同士で反転していない時と反転している時のスコア
map<pair<size_t,size_t>,pair<int,int>> VCFCollection::compute_scores() const {
	// VCFが変わる部分を抽出する
	const vector<Joint>	joints = extract_joints();
	
	map<pair<size_t,size_t>,pair<int,int>>	dic_scores;
	for(auto p = joints.begin(); p != joints.end(); ++p) {
		const auto&	pair1 = p->first;
		const auto&	pair2 = p->second;
		const size_t&	i1 = pair1.first;
		const size_t&	i2 = pair2.first;
		pair<size_t,size_t>	key(i1, i2);	// pair of vcf indices
		dic_scores[key].first  += score_joint(*p, true);
		dic_scores[key].second += score_joint(*p, false);
	}
	return dic_scores;
}

// ヘテロ親のハプロタイプを逆にしたときとしていないときの繋がり具合
// 各サンプルでヘテロ親のどちらから来たかが変わっていれば+1
int VCFCollection::score_connection(const vector<bool>& bs,
					const map<VCF_INDEX_PAIR,pair<int,int>>& scores) const {
	int score = 0;
	for(auto p = scores.begin(); p != scores.end(); ++p) {
		const size_t	i1 = p->first.first;
		const size_t	i2 = p->first.second;
		const auto&		v = p->second;
		score += bs[i1] == bs[i2] ? v.first : v.second;
	}
	return score;
}

void VCFCollection::determine_haplotype() {
	if(vcfs.size() == 1U)
		return;
	
	// vcfsが最も矛盾なくつながるように親を反転する
	const map<pair<size_t,size_t>,pair<int,int>>	scores = compute_scores();
	vector<vector<bool>>	bss = make_boolean_vectors(vcfs.size());
	vector<bool>	min_bs = bss.front();
	int	min_score = score_connection(min_bs, scores);
	for(auto p = bss.begin() + 1; p != bss.end(); ++p) {
		const int	score = score_connection(*p, scores);
		if(score < min_score) {
			min_bs = *p;
			min_score = score;
		}
	}
	
	for(size_t i = 0; i < vcfs.size(); ++i) {
		if(min_bs[i])
			vcfs[i]->inverse_haplotype();
	}
	
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		(*p)->update_genotypes();
}

VCFFamily *VCFCollection::join() const {
	vector<VCFFamilyRecord *>	records(max_index() + 1, NULL);
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*subvcf = *p;
		for(size_t i = 0U; i < subvcf->size(); ++i) {
			VCFImputableRecord	*record = subvcf->get_record(i);
			records[record->get_index()] = record->copy();
		}
	}
	return new VCFFamily(vcfs.front()->get_header(),
							vcfs.front()->get_samples(), records);
}

VCFFamily *VCFCollection::impute_one_parent(VCFHeteroHomo *vcf,
										const vector<const Map *>& chr_maps,
										bool is_mat, int MIN_CROSSOVER, int T,
										BiasProbability *bias_probability) {
	OptionImpute	op(4, 20, 5);
	vector<VCFFamily *>	vcfs;
	vector<VCFHeteroHomo *>	chr_vcfs = vcf->divide_into_chromosomes(chr_maps);
	for(size_t i = 0U; i < chr_vcfs.size(); ++i) {
		VCFHeteroHomo	*subvcf1 = chr_vcfs[i];
		auto	*collection = VCFImputable::determine_haplotype(subvcf1,
												&op, is_mat, bias_probability);
		delete subvcf1;
		collection->impute(op.min_positions, MIN_CROSSOVER, T);
		if(collection->empty()) {
			delete collection;
			continue;
		}
		collection->determine_haplotype();
		auto	*subvcf2 = collection->join();
		delete collection;
		vcfs.push_back(subvcf2);
	}
	auto	*new_vcf = VCFFamily::join(vcfs);
	Common::delete_all(vcfs);
	vcf->copy_chrs(new_vcf);
	return new_vcf;
}

VCFFamily *VCFCollection::impute_family_vcf(VCFHeteroHomo *mat_vcf,
											VCFHeteroHomo *pat_vcf,
											const vector<const Map *>& chr_maps,
											int MIN_CROSSOVER, int T) {
	auto	*bias_probability = new BiasProbability(0.01);
	VCFFamily	*impute_mat_vcf = impute_one_parent(mat_vcf, chr_maps,
													true, MIN_CROSSOVER, T,
													bias_probability);
	VCFFamily	*impute_pat_vcf = impute_one_parent(pat_vcf, chr_maps,
													false, MIN_CROSSOVER, T,
													bias_probability);
	VCFFamily	*vcf =  VCFFamily::merge(impute_mat_vcf, impute_pat_vcf);
	delete impute_mat_vcf;
	delete impute_pat_vcf;
	delete bias_probability;
	return vcf;
}
