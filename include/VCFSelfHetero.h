#ifndef __VCFSELFHETERO
#define __VCFSELFHETERO

#include "VCFGeno.h"
#include "VCFSelfHeteroRecord.h"
#include "Map.h"
#include "invgraph.h"

class Map;
class PedigreeTable;
class Family;
class VCFOriginal;
class VCFSelfHetero;
class VCFSelfHeteroRecord;
class InverseGraph;
class Option;
class OptionImpute;

typedef std::pair<std::vector<VCFSelfHetero *>,
						std::vector<VCFSelfHeteroRecord *>>	ImpResult;


//////////////////// VCFSelfHetero ////////////////////

class VCFSelfHetero : public VCFGenoBase, public VCFMeasurable {
protected:
	std::vector<VCFSelfHeteroRecord *>	records;
	
public:
	VCFSelfHetero(const STRVEC& s,
					const std::vector<VCFSelfHeteroRecord *>& rs,
					const Map& m, const VCFSmall *vcf);
	~VCFSelfHetero();
	
	///// virtual methods for VCFGenoBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// non-virtual methods /////
	std::size_t num_progenies() const { return num_samples() - 1; }
	const std::vector<VCFSelfHeteroRecord *>& get_records() const {
		return records;
	}
	
	void set_records(const std::vector<VCFSelfHeteroRecord *>& rs) {
		records = rs;
	}
	void clear_records() { records.clear(); }
	
	// decide heterozygous parent haplotypes
	std::vector<int> make_parent_haplotypes(const InvGraph& graph) const;
	std::vector<int> create_haplotype(std::size_t v0,
										const InvGraph& graph) const;
	VCFSelfHetero *make_subvcf(const InvGraph& graph) const;
	std::pair<std::vector<VCFSelfHetero *>, std::vector<VCFSelfHeteroRecord *>>
						determine_haplotype(const OptionImpute *option) const;
	std::pair<int,bool> distance(const std::vector<int>& gts1,
							const std::vector<int>& gts2, int max_dist) const;
	std::pair<std::vector<VCFSelfHetero *>, std::vector<VCFSelfHeteroRecord *>>
														impute(int num_threads);
	
private:
	double record_cM(std::size_t i) const { return cM(records[i]->get_pos()); }
	InvGraph make_graph(double max_dist) const;
	std::string make_seq(std::size_t i) const;
	std::string impute_each_sample_seq(std::size_t i,
								const std::vector<double>& cMs, double min_c);
	void impute_each_sample(std::size_t i,
								const std::vector<double>& cMs, double min_c);
	void impute_each(const OptionImpute *option);
	const OptionImpute *create_option(int num_threads) const;
	
private:
	static std::pair<double, bool> distance(const std::vector<int>& gts1,
											const std::vector<int>& gts2);
	static double dist_with_NA(int right, int counter_NA, int N);
	static void impute_in_thread(void *config);
};
#endif
