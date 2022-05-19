#ifndef __VCFIMPUTABLE
#define __VCFIMPUTABLE

#include "VCFFamily.h"
#include "VCFHeteroHomo.h"
#include "graph.h"
#include "Baum_Welch_with_fixed_Ts.h"

class Map;
class BiasProbability;
class PedigreeTable;
class VCFOriginal;
class VCFCollection;
class OptionImpute;


//////////////////// VCFImputableRecord ////////////////////

class VCFImputableRecord : public VCFHeteroHomoRecord {
	int		haplo[2];
	bool	mat_hetero;
	std::vector<int>	int_gts;
	std::size_t	index;	// index in VCF
	int	group_id;		// index of graph
	std::vector<int>	which_comes_from;
	
public:
	VCFImputableRecord(const STRVEC& v, const STRVEC& s, int h1, int h2,
						bool b, std::size_t i, const std::vector<int>& gts);
	
	std::size_t get_index() const { return index; }
	int mat_int_gt() const { return int_gts[0]; }
	int pat_int_gt() const { return int_gts[1]; }
	STRVEC get_GTs() const;
	std::string homo() const { return samples[mat_hetero?1:0]; }
	bool is_mat_hetero() const { return mat_hetero; }
	int homo_int_gt() const { return mat_hetero ? pat_int_gt() : mat_int_gt(); }
	int get_which_comes_from(size_t i) const { return which_comes_from[i]; }
	int difference(const VCFImputableRecord *other, bool inv) const;
	
	void set_index(size_t i) { index = i; }
	void set_group_id(int id) { group_id = id; }
	void which_comes_from_hetero_parent_chromosomes();
	
	void set_int_gts_by_which_comes_from(const std::vector<int>& hs);
	void update_genotypes();
	void inverse_haplotype();
	
private:
	int which_comes_from_hetero_parent_chromosomes_each(int gt);
	int recover_int_gt(int h);
	int GT_position() const;
	
	std::string decode_gt(int h);
	std::string decode_parent_gt(int int_gt);
	std::string update_gt(const std::string& gt,
							const std::string& s, int gt_pos);

};


//////////////////// VCFImputable ////////////////////

class VCFImputable : public VCFHeteroHomo {
	typedef std::map<std::pair<char,char>,double>	Matrix;
	typedef std::vector<std::map<char,double>>	Table;
	
	struct ConfigThread {
		const VCFImputable	*vcf;
		const std::vector<Matrix>&	Ts;
		const std::size_t	first;
		const double	MIN;
		const int	num_thread;
		std::vector<std::string>& hs;
		
		ConfigThread(const VCFImputable* v, const std::vector<Matrix>& mat,
						int f, double m, int n, std::vector<std::string>& h) :
					vcf(v), Ts(mat), first(f), MIN(m), num_thread(n), hs(h) { }
	};
	
	std::vector<VCFImputableRecord *>	imp_records;
	const bool	mat_hetero;
	
public:
	VCFImputable(const std::vector<STRVEC>& h, const STRVEC& s,
						std::vector<VCFImputableRecord *> rs, const Map& m);
	~VCFImputable() { }
	
	VCFImputableRecord *get_record(size_t i) const { return imp_records[i]; }
	std::vector<STRVEC> get_GT_table() const;
	std::size_t max_index() const;
	
	void set_group_id(int id);
	void impute(double MIN_CROSSOVER, int T);
	void update_genotypes();
	void inverse_haplotype();
	VCFCollection *divide(const OptionImpute *option, bool is_mat) const;
	
private:
	std::vector<VCFHeteroHomoRecord *> to_VCFHeteroHomoRecord(
										std::vector<VCFImputableRecord *>& rs);
	Graph::WeightedGraph trim_inverse(const Graph::InvGraph& graph) const;
	Graph::InvGraph filter_graph(const Graph::InvGraph& graph,
									const Graph::WeightedGraph& tree) const;
	bool is_all_same_without_N(const std::string& seq) const;
	std::string create_same_color_string(const std::string& seq) const;
	std::string impute_each_seq(const std::string& seq,
								const std::vector<Matrix>& Ts,
								double MIN_CROSSOVER) const;
	std::vector<Matrix> create_transition_probability_matrix() const;
	std::string make_seq(std::size_t i) const;
	void which_comes_from_hetero_parent_chromosomes();
	
public:
	static void walk(std::size_t v0, std::vector<std::vector<int>>& haplo,
										const Graph::InvGraph& tree_graph);
	static Graph::InvGraph make_graph(VCFHeteroHomo *vcf, int max_dist);
	static VCFImputable *make_subvcf(VCFHeteroHomo *vcf,
									const Graph::InvGraph& graph, bool is_mat,
									BiasProbability *bias_probability);
	static void renumber_indices(std::vector<VCFImputable *>& vcfs);
	static VCFCollection *determine_haplotype(VCFHeteroHomo *vcf,
									const OptionImpute *option, bool is_mat,
									BiasProbability *bias_probability);
	static VCFImputable *create(const std::vector<STRVEC>& header,
								const std::vector<VCFHeteroHomoRecord *>& rs,
								bool hatero,
								const std::vector<std::vector<int>>& haplo,
								const std::vector<std::size_t>& indices,
								const Map& m);
	static bool determine_which_parents_is_hetero(VCFFamily *vcf);
	static std::vector<std::vector<int>> make_parent_haplotypes(std::size_t L,
												const Graph::InvGraph& graph);
	static void impute_by_thread(void *config);
};

#endif
