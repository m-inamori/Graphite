#ifndef __VCFCOLLECTION
#define __VCFCOLLECTION

#include "VCFImputable.h"

using namespace std;

class Map;
class BiasProbability;


//////////////////// VCFCollection ////////////////////

class VCFCollection {
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
	std::vector<std::pair<std::size_t,std::size_t>>
	make_which_vcf_table() const;
	int score_connection(const std::vector<bool>& bs) const;
	
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
