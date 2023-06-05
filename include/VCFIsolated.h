#ifndef __VCFISOLATED
#define __VCFISOLATED

#include "VCF.h"
#include "Map.h"


//////////////////// VCFIsolated ////////////////////

// imputed samples at the beginning, followed by reference samples
class VCFIsolated : public VCFSmall, VCFMeasurable {
public:
	typedef std::pair<std::size_t, std::size_t>	Haplotype;
	typedef std::pair<Haplotype, Haplotype>	HaplotypePair;
	
	struct ConfigThread {
		const std::vector<VCFIsolated *>&	vcfs;
		const std::size_t	first;
		const int	num_threads;
		
		ConfigThread(const std::vector<VCFIsolated *>& vcfs_,
											std::size_t f, int n) :
									vcfs(vcfs_), first(f), num_threads(n) { }
		
		std::size_t size() const { return vcfs.size(); }
	};
	
protected:
	const std::size_t	num_imputed_samples;
	
public:
	VCFIsolated(const std::vector<STRVEC>& h, const STRVEC& s,
				std::vector<VCFRecord *> rs, std::size_t nis, const Map& m);
	~VCFIsolated() { }
	
	void impute();
	VCFSmall *extract_isolated_samples() const;
	
private:
	bool is_block(const VCFRecord *record,
					const std::vector<VCFRecord *>& rs) const;
	
	// divide VCF by 1 cM
	std::vector<VCFIsolated *> divide_by_cM() const;
	
	int get_single_gt(const VCFRecord *record, Haplotype hap) const;
	int score_each(Haplotype hap1, Haplotype hap2,
							std::size_t i, VCFRecord *record) const;
	int score(Haplotype hap1, Haplotype hap2, std::size_t i) const;
	std::vector<std::pair<Haplotype, Haplotype>>
				collect_optimal_haplotype_pairs(std::size_t i) const;
	void set_haplotype(HaplotypePair hap, std::size_t i);
	int match_score(HaplotypePair prev_hap, HaplotypePair hap) const;
	std::vector<HaplotypePair> collect_max_score(
									const std::vector<HaplotypePair>& combs,
									HaplotypePair prev_hap) const;
	HaplotypePair impute_cM_each_sample(HaplotypePair prev_hap, std::size_t i);
	std::vector<HaplotypePair> impute_cM(
								const std::vector<HaplotypePair>& prev_haps);
	
public:
	static std::vector<VCFIsolated *>
			create(const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
					const std::vector<std::string>& samples,
					const std::vector<std::string>& references,
					const Map& gmap, int num_threads);
	static std::vector<VCFSmall *> impute_all(
										const std::vector<VCFIsolated *>& vcfs,
										int num_threads);
	
private:
	static std::vector<std::vector<std::size_t>>
			divide_columns(const std::vector<std::size_t>& cs, int num);
	static void impute_parellel(const std::vector<VCFIsolated *>& vcfs,
															int num_threads);
};
#endif
