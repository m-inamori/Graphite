#ifndef __VCFHETEROHOMO
#define __VCFHETEROHOMO

#include "VCFImpFamily.h"
#include "option.h"
#include "graph.h"
#include "Baum_Welch_with_fixed_Ts.h"

class Map;
class BiasProbability;
class PedigreeTable;
class VCFOriginal;
class VCFHeteroHomoRecord;
class VCFHeteroHomo;

typedef std::pair<std::vector<VCFHeteroHomo *>,
						std::vector<VCFHeteroHomoRecord *>>	ImpResult;


//////////////////// VCFHeteroHomoRecord ////////////////////

class VCFHeteroHomoRecord : public VCFImpFamilyRecord {
	std::vector<int>	which_comes_from;
	
public:
	VCFHeteroHomoRecord(const STRVEC& v, const STRVEC& s,
							int i, WrongType type, ParentComb c) :
							VCFImpFamilyRecord(v, s, i, type, c),
							which_comes_from(num_progenies(), -1) { }
	
	int get_which_comes_from(int i) const { return which_comes_from[i]; }
	bool is_mat_hetero() const { return this->mat_int_gt() == 1; }
	bool is_pat_hetero() const { return this->pat_int_gt() == 1; }
	
	bool is_homohomo() const { return false; }
	bool is_imputable() const {
		return this->wrong_type == WrongType::RIGHT;
	}
	FillType get_fill_type() const;
	
	std::vector<int> genotypes_from_hetero_parent() const;
	void set_haplo(int h);
	void set_int_gt_by_which_comes_from(const std::vector<int>& ws);
};


//////////////////// VCFHeteroHomo ////////////////////

class VCFHeteroHomo : public VCFFamily {
public:
	typedef std::pair<std::string,std::string>	Parents;
	
protected:
	std::vector<VCFHeteroHomoRecord *>	hh_records;
	const Map&	genetic_map;
	
public:
	VCFHeteroHomo(const std::vector<STRVEC>& h, const STRVEC& s,
						std::vector<VCFHeteroHomoRecord *> rs, const Map& m);
	~VCFHeteroHomo() { }
	
	bool is_mat_hetero() const { return hh_records.front()->is_mat_hetero(); }
	
	double cM(std::size_t i) const;
	void set_records(const std::vector<VCFHeteroHomoRecord *>& rs);
	void set_records_base(const std::vector<VCFHeteroHomoRecord *>& rs);
	void clear_records();
	
	// haplotype1の各レコードのGenotypeが0なのか1なのか
	std::vector<int> make_parent_haplotypes(const Graph::InvGraph& graph) const;
	VCFHeteroHomo *make_subvcf(const Graph::InvGraph& graph) const;
	std::pair<std::vector<VCFHeteroHomo *>, std::vector<VCFHeteroHomoRecord *>>
						determine_haplotype(const OptionImpute *option) const;
	std::pair<std::vector<VCFHeteroHomo *>, std::vector<VCFHeteroHomoRecord *>>
																	impute();
	// 共通のヘテロ親はどれだけマッチしているか
	std::pair<int, int> match(const VCFHeteroHomo *other) const;
	void inverse_hetero_parent_phases();
	
	VCFHeteroHomoRecord *get_record(std::size_t i) const {
		return hh_records[i];
	}
	const std::vector<VCFHeteroHomoRecord *>& get_records() const {
		return hh_records;
	}
	const Map& get_map() const { return genetic_map; }
	std::pair<int,bool> distance(const std::vector<int>& gts1,
							const std::vector<int>& gts2, int max_dist) const;
	
private:
	Graph::InvGraph make_graph(double max_dist) const;
	std::string make_seq(std::size_t i) const;
	std::string impute_each_sample_seq(int i,
								const std::vector<double>& cMs, double min_c);
	void impute_each(const OptionImpute *option);
	const OptionImpute *create_option() const;
	
public:
	// ヘテロ親が同じVCFを集めて補完する
	// ついでにphaseもなるべく同じになるように変更する
	static ImpResult impute_vcfs(const std::vector<VCFHeteroHomo *>& vcfs,
											const Option* op, int num_threads);
	static void inverse_phases(const std::vector<VCFHeteroHomo *>& vcfs);
	static std::vector<std::vector<std::tuple<int, int, int>>> make_vcf_graph(
										const std::vector<VCFHeteroHomo *>& vcfs);
	
private:
	static std::pair<double, bool> distance(const std::vector<int>& gts1,
											const std::vector<int>& gts2);
	static double dist_with_NA(int right, int counter_NA, int N);
	static void walk(std::size_t v0, std::vector<int>& haplo,
									const Graph::InvGraph& tree_graph);
	static bool is_all_same_without_N(const std::string& seq);
	static std::string create_same_color_string(const std::string& seq);
	static std::vector<std::vector<bool>> enumerate_bools(std::size_t L);
	static int match_score(const std::vector<bool>& invs,
			const std::vector<std::vector<std::tuple<int, int, int>>>& graph);
	static std::vector<bool> optimize_phase_inversions(
			const std::vector<std::vector<std::tuple<int, int, int>>>& graph);
};
#endif
