#ifndef __VCFCOLLECTION
#define __VCFCOLLECTION

#include "VCFImputable.h"

using namespace std;

class Map;
class BiasProbability;


//////////////////// VCFCollection ////////////////////

class VCFCollection {
	typedef std::pair<std::pair<std::size_t, std::size_t>,
					  std::pair<std::size_t, std::size_t>>	Joint;
	typedef std::pair<std::size_t,std::size_t>	VCF_INDEX_PAIR;
	
private:
	std::vector<VCFImputable *>	vcfs;
	
public:
	VCFCollection(const std::vector<VCFImputable *>& vcfs_);
	~VCFCollection();
	
	bool empty() const { return vcfs.empty(); }
	std::size_t max_index() const;
	
	void impute(int min_positions, double MIN_CROSSOVER, int T);
	void determine_haplotype();
	VCFFamily *join() const;
	
private:
	std::vector<std::vector<bool>> all_boolean_vectors(std::size_t n) const;
	std::vector<std::vector<bool>> make_random_boolean_vectors(
											size_t num_vec, size_t n) const;
	std::vector<VCF_INDEX_PAIR> make_which_vcf_table() const;
	int score_connection(const std::vector<bool>& bs,
			const std::map<VCF_INDEX_PAIR,std::pair<int,int>>& scores) const;
	std::vector<Joint> extract_joints() const;
	std::map<std::pair<std::size_t,std::size_t>,std::pair<int,int>>
	compute_scores() const;
	int score_joint(const Joint& joint, bool b) const;
	
public:
	static VCFFamily *impute_one_parent(VCFHeteroHomo *vcf,
									const std::vector<const Map *>& chr_maps,
									bool is_mat, int MIN_CROSSOVER, int T,
									BiasProbability *bias_probability);
	static VCFFamily *impute_family_vcf(VCFHeteroHomo *mat_vcf,
										VCFHeteroHomo *pat_vcf,
										const vector<const Map *>& chr_maps,
										int MIN_CROSSOVER, int T);
};
#endif
