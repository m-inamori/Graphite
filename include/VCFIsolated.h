#ifndef __VCFISOLATED
#define __VCFISOLATED

#include "VCF.h"
#include "Map.h"
#include "VCFClippable.h"

class OptionSmall;


//////////////////// VCFIsolated ////////////////////

// imputed samples at the beginning, followed by reference samples
class VCFIsolated : public VCFClippable {
public:
	struct ConfigThread {
		const std::vector<VCFIsolated *>&	vcfs;
		const std::size_t	first;
		const int	num_threads;
		
		ConfigThread(const std::vector<VCFIsolated *>& vcfs_,
											std::size_t f, int n) :
									vcfs(vcfs_), first(f), num_threads(n) { }
		
		std::size_t size() const { return vcfs.size(); }
	};
	
private:
	std::vector<GenoRecord *>	records;
	const std::size_t	num_imputed_samples;
	
public:
	VCFIsolated(const STRVEC& s, std::size_t nis,
				std::vector<GenoRecord *>& rs,
				const Map& m, const VCFSmall *vcf);
	~VCFIsolated();
	
	///// virtual methods for VFSmallBase /////
	std::size_t size() const override { return records.size(); }
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFClippable /////
	std::vector<Haplotype> collect_haplotypes_mat(
									std::size_t sample_index) const override;
	std::vector<Haplotype> collect_haplotypes_pat(
									std::size_t sample_index) const override;
	void set_gts(const std::vector<std::string>& gts, std::size_t sample_index);
	
	///// non-virtual methods /////
	void impute();
	void add_record(GenoRecord *record) { records.push_back(record); }
	
	std::vector<Haplotype> collect_haplotype_from_refs() const;
	
private:
	std::vector<VCFIsolated *> divide_by_cM() const;
	std::vector<HaplotypePair> impute_cM(
								const std::vector<HaplotypePair>& prev_haps);
	
public:
	static VCFIsolated *create(const VCFSmall *orig_vcf,
								const VCFGeno *merged_vcf,
								const STRVEC& samples,
								const STRVEC& references,
								const OptionSmall& op);
};
#endif
