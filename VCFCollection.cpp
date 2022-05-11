#include <random>
#include <cassert>
#include "VCFCollection.h"
#include "option.h"
#include "Map.h"

using namespace std;


//////////////////// VCFCollection ////////////////////

VCFCollection::VCFCollection(const vector<VCFImputable *>& vcfs_) :
														vcfs(vcfs_) {
	for(int i = 0; i < (int)vcfs.size(); ++i)
		vcfs[i]->set_group_id(i);
}

VCFCollection::~VCFCollection() {
	// VCFHeteroHomo::divide_into_chromosomesで作ったMapはvcfsで共用
	delete &(vcfs.front()->get_map());
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

vector<vector<bool>> VCFCollection::all_boolean_vectors(size_t n) const {
	// D&C
	if(n == 1U) {
		vector<vector<bool>>	bs = { vector<bool>(1, true),
										vector<bool>(1, false) };
		return bs;
	}
	else if(n % 2 == 1U) {
		const auto	part_vectors = all_boolean_vectors(n - 1);
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
		const auto	part_vectors = all_boolean_vectors(n / 2);
		vector<vector<bool>>	vectors;
		for(auto p = part_vectors.begin(); p != part_vectors.end(); ++p) {
			vector<bool>	v1 = *p;
			for(auto q = part_vectors.begin(); q != part_vectors.end(); ++q) {
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
											size_t num_vec, size_t n) const {
	std::mt19937	mt(123456789);
	std::uniform_int_distribution<int>	dist(0, 1);
	
	vector<vector<bool>>	vectors(num_vec, vector<bool>(n));
	for(size_t i = 0; i < num_vec; ++i) {
		for(size_t j = 0; j < n; ++j)
			vectors[i][j] = dist(mt) == 0;
	}
	return vectors;
}

vector<pair<size_t,size_t>> VCFCollection::make_which_vcf_table() const {
	const size_t	num_markers = max_index() + 1;
	vector<pair<size_t,size_t>>	which_vcfs(num_markers);
	for(size_t i = 0U; i < vcfs.size(); ++i) {
		const auto	*vcf = vcfs[i];
		for(size_t j = 0U; j < vcfs.size(); ++j) {
			const auto	*record = vcf->get_record(j);
			which_vcfs[record->get_index()] = pair<size_t,size_t>(i, j);
		}
	}
	return which_vcfs;
}

// ヘテロ親のハプロタイプを逆にしたときとしていないときの繋がり具合
// 各サンプルでヘテロ親のどちらから来たかが変わっていれば+1
int VCFCollection::score_connection(const vector<bool>& bs) const {
	const auto	which_vcfs = make_which_vcf_table();
	int score = 0;
	for(auto p = which_vcfs.begin(); p != which_vcfs.end() - 1; ++p) {
		const size_t	i1 = p->first;
		const size_t	j1 = p->second;
		const size_t	i2 = (p+1)->first;
		const size_t	j2 = (p+2)->second;
		if(i1 != i2) {
			const auto	*record1 = vcfs[i1]->get_record(j1);
			const auto	*record2 = vcfs[i2]->get_record(j2);
			score += record1->difference(record2, bs[i1] ^ bs[i2]);
		}
	}
	return score;
}

// あとでDPにする
// そうすれば乱数を使わなくて済む
void VCFCollection::determine_haplotype() {
	// vcfsが最も矛盾なくつながるように親を反転する
	const size_t	bs_limit = 10;
	vector<vector<bool>>	bss;
	if(vcfs.size() <= bs_limit) {
		bss = all_boolean_vectors(vcfs.size());
	}
	else {
		bss = make_random_boolean_vectors(1U << bs_limit, vcfs.size());
	}
	
	vector<bool>	min_bs = bss.front();
	int	min_score = score_connection(min_bs);
	for(auto p = bss.begin() + 1; p != bss.end(); ++p) {
		const int	score = score_connection(*p);
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

VCFFamily *VCFCollection::impute_one_parent(VCFHeteroHomo *vcf, const Map& gmap,
										bool is_mat, int MIN_CROSSOVER, int T) {
	OptionImpute	op(4, 20, 5);
	vector<VCFFamily *>	vcfs;
	vector<VCFHeteroHomo *>	chr_vcfs = vcf->divide_into_chromosomes();
	for(size_t i = 0U; i < chr_vcfs.size(); ++i) {
		VCFHeteroHomo	*subvcf1 = chr_vcfs[i];
		auto	*collection = VCFImputable::determine_haplotype(subvcf1,
																&op, is_mat);
		delete subvcf1;
		collection->impute(op.min_positions, MIN_CROSSOVER, T);
		if(collection->empty())
			continue;
		auto	*subvcf2 = collection->join();
		delete collection;
		vcfs.push_back(subvcf2);
	}
	auto	*new_vcf = VCFFamily::join(vcfs);
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		delete *p;
	vcf->copy_chrs(new_vcf);
	return new_vcf;
}

VCFFamily *VCFCollection::impute_family_vcf(VCFHeteroHomo *mat_vcf,
											VCFHeteroHomo *pat_vcf,
											const Map& gmap,
											int MIN_CROSSOVER, int T) {
	VCFFamily	*impute_mat_vcf = impute_one_parent(mat_vcf, gmap,
													true, MIN_CROSSOVER, T);
	VCFFamily	*impute_pat_vcf = impute_one_parent(pat_vcf, gmap,
													false, MIN_CROSSOVER, T);
	return VCFFamily::merge(impute_mat_vcf, impute_pat_vcf);
}
