#ifndef __VCFONEPARENTPHASED
#define __VCFONEPARENTPHASED

#include "VCFFamily.h"
#include "Map.h"
#include "option.h"
#include "graph.h"

class PedigreeTable;
class Family;
class VCFOriginal;


//////////////////// VCFOneParentPhased ////////////////////

class VCFOneParentPhased : public VCFFamily, VCFMeasurable {
public:
	struct ConfigThread {
		const VCFSmall	*orig_vcf;
		const VCFSmall	*merged_vcf;
		const std::vector<std::pair<const Family *, bool>>&	families;
		const Map& geno_map;
		const std::size_t	first;
		const int	num_threads;
		std::vector<VCFFamily *>&	results;
		
		ConfigThread(const VCFSmall *o, const VCFSmall *m,
					 const std::vector<std::pair<const Family *, bool>>& fs,
					 const Map& gmap, std::size_t f,
					 int T, std::vector<VCFFamily *>& rs) :
			 						orig_vcf(o), merged_vcf(m),
			 						families(fs), geno_map(gmap),
			 						first(f), num_threads(T), results(rs) { }
		
		std::size_t size() const { return families.size(); }
	};
	
private:
	const bool	is_mat_phased;
	
public:
	VCFOneParentPhased(const std::vector<STRVEC>& h, const STRVEC& s,
										std::vector<VCFFamilyRecord *> rs,
										bool is_mat_phased, const Map& m);
	~VCFOneParentPhased() { }
	
	bool is_mat_hetero() const;
	
	void impute();
	
private:
	char determine_which_comes_from(VCFFamilyRecord *record,
												std::size_t i) const;
	double record_cM(std::size_t i) const { return cM(records[i]->pos()); }
	std::string make_seq(std::size_t i) const;
	std::string impute_sample_seq(std::size_t i,
							const std::vector<double>& cMs, double min_c) const;
	void update_each(std::size_t i, std::size_t j, char c);
	void update(std::size_t i, const std::string& seq);
	static void impute_in_thread(void *config);
	
public:
	static VCFFamily *impute_by_parent(const VCFSmall *orig_vcf,
										const VCFSmall *parent_imputed_vcf,
										const STRVEC& samples,
										bool is_mat_phased, const Map& gmap);
	static std::vector<VCFFamily *> impute_all_by_parent(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				const std::vector<std::pair<const Family *, bool>>& families,
				const Map& geno_map, int num_threads);
	// samplesに対応するRecordを作る
	static VCFRecord *merge_records(const std::vector<VCFFamily *>& vcfs,
										std::size_t i, const STRVEC& samples);
};
#endif
