#ifndef __VCFHETEROHOMO
#define __VCFHETEROHOMO

#include "VCFImpFamilyRecord.h"
#include "Map.h"
#include "invgraph.h"
#include "Baum_Welch_with_fixed_Ts.h"

class Map;
class BiasProbability;
class PedigreeTable;
class Family;
class VCFOriginal;
class VCFHeteroHomoRecord;
class VCFHeteroHomo;
class InverseGraph;
class Option;
class OptionImpute;

typedef std::pair<std::vector<VCFHeteroHomo *>,
						std::vector<VCFHeteroHomoRecord *>>	ImpResult;


//////////////////// VCFHeteroHomoRecord ////////////////////

class VCFHeteroHomoRecord : public VCFImpFamilyRecord {
	std::vector<int>	which_comes_from;
	
public:
	VCFHeteroHomoRecord(ll pos, const std::vector<int>& geno,
							int i, WrongType type, ParentComb c) :
							VCFImpFamilyRecord(pos, geno, i, type, c),
							which_comes_from(num_progenies(), -1) { }
	
	int get_which_comes_from(int i) const { return which_comes_from[i]; }
	
	bool is_homohomo() const override { return false; }
	bool is_imputable() const override {
		return this->wrong_type == WrongType::RIGHT;
	}
	FillType get_fill_type() const override;
	
	std::vector<int> alleles_from_hetero_parent() const;
	void set_haplo(int h);
	void set_int_gt_by_which_comes_from(int w, std::size_t i);
};


//////////////////// VCFHeteroHomo ////////////////////

class VCFHeteroHomo : public VCFFamilyBase, public VCFMeasurable {
public:
	struct ConfigThreadCleanSeq {
		const std::size_t	first;
		const std::size_t	num_threads;
		const std::vector<double>&	cMs;
		const double	min_c;
		VCFHeteroHomo	*vcf;
		
		ConfigThreadCleanSeq(std::size_t i, std::size_t	n,
							 const std::vector<double>& cMs_,
							 double m, VCFHeteroHomo *vcf_) :
					first(i), num_threads(n), cMs(cMs_), min_c(m), vcf(vcf_) { }
	};
	
	typedef std::pair<std::string,std::string>	Parents;
	
protected:
	std::vector<VCFHeteroHomoRecord *>	records;
	
public:
	VCFHeteroHomo(const STRVEC& s,
					const std::vector<VCFHeteroHomoRecord *>& rs,
					const Map& m, const VCFSmall *vcf);
	~VCFHeteroHomo();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	const std::vector<VCFHeteroHomoRecord *>& get_records() const {
		return records;
	}
	bool is_mat_hetero() const { return records.front()->is_mat_hetero(); }
	
	void set_records(const std::vector<VCFHeteroHomoRecord *>& rs) {
		records = rs;
	}
	void clear_records() { records.clear(); }
	
	// decide heterozygous parent haplotypes
	std::vector<int> make_parent_haplotypes(const InvGraph& graph) const;
	std::vector<int> create_haplotype(std::size_t v0,
										const InvGraph& graph) const;
	VCFHeteroHomo *make_subvcf(const InvGraph& graph) const;
	std::pair<std::vector<VCFHeteroHomo *>, std::vector<VCFHeteroHomoRecord *>>
						determine_haplotype(const OptionImpute *option) const;
	std::pair<std::vector<VCFHeteroHomo *>, std::vector<VCFHeteroHomoRecord *>>
														clean(int num_threads);
	// How well do the common heteroparents match?
	std::pair<int, int> match(const VCFHeteroHomo *other) const;
	void inverse_hetero_parent_phases();
	
	std::pair<int,bool> distance(const std::vector<int>& gts1,
							const std::vector<int>& gts2, int max_dist) const;
	
private:
	double record_cM(std::size_t i) const { return cM(records[i]->get_pos()); }
	InvGraph make_graph(double max_dist) const;
	std::string make_seq(std::size_t i) const;
	std::string clean_each_sample_seq(std::size_t i,
								const std::vector<double>& cMs, double min_c);
	void clean_each_sample(std::size_t i,
								const std::vector<double>& cMs, double min_c);
	void clean_each_vcf(const OptionImpute *option);
	const OptionImpute *create_option(int num_threads) const;
	
public:
	// Create a pair of VCFHeteroHomo for each Family
	// and store it for each parent
	static std::tuple<VCFHeteroHomo *, VCFHeteroHomo *,
						std::vector<VCFHeteroHomoRecord *>>
		make_VCFHeteroHomo(const std::vector<VCFHeteroHomoRecord *>& records,
						   const STRVEC& samples, const Map& geno_map,
						   const VCFSmall *vcf);
	// collect VCFs whose hetero parent is the same, and impute it
	// and change the phase to be as same as possible
	static ImpResult clean_vcfs(
						const std::vector<VCFHeteroHomoRecord *>& records,
						const STRVEC& samples, const Map& geno_map,
						int num_threads, const VCFSmall *vcf);
	static void inverse_phases(const std::vector<VCFHeteroHomo *>& vcfs);
	static const InverseGraph *make_vcf_graph(
										const std::vector<VCFHeteroHomo *>& vcfs);
	
private:
	static std::pair<double, bool> distance(const std::vector<int>& gts1,
											const std::vector<int>& gts2);
	static double dist_with_NA(int right, int counter_NA, int N);
	static std::vector<bool> optimize_phase_inversions(
			const std::vector<std::vector<std::tuple<int, int, int>>>& graph);
	static void clean_in_thread(void *config);
};
#endif
