#ifndef __VCFONEPARENTPHASED
#define __VCFONEPARENTPHASED

#include "VCFFamily.h"
#include "Map.h"
#include "VCFImputable.h"
#include "option.h"
#include "graph.h"

class PedigreeTable;
class Family;
class VCFOriginal;


//////////////////// VCFOneParentPhased ////////////////////

class VCFOneParentPhased : public VCFBase, public VCFFamilyBase,
											public VCFImputable {
public:
	struct ConfigThread {
		const std::vector<VCFOneParentPhased *>	vcfs;
		const std::size_t	first;
		const int	num_threads;
		
		ConfigThread(const std::vector<VCFOneParentPhased *>& vcfs_,
												std::size_t f, int T) :
			 						vcfs(vcfs_), first(f), num_threads(T) { }
		
		std::size_t size() const { return vcfs.size(); }
	};
	
private:
	std::vector<VCFFamilyRecord *>	records;
	const bool		is_mat_phased;
	const VCFSmall	*ref_vcf;
	
public:
	VCFOneParentPhased(const std::vector<STRVEC>& h, const STRVEC& s,
									const std::vector<VCFFamilyRecord *>& rs,
									bool is_mat_phased, const Map& m,
									const VCFSmall *ref);
	VCFOneParentPhased(const VCFOneParentPhased&) = delete;
	VCFOneParentPhased& operator=(const VCFOneParentPhased&) = delete;
	~VCFOneParentPhased();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	VCFFamilyRecord *get_family_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFImputable /////
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
	std::vector<Haplotype> collect_haplotypes_mat(
									std::size_t sample_index) const override;
	std::vector<Haplotype> collect_haplotypes_pat(
									std::size_t sample_index) const override;
	
	///// non-virtual methods /////
	bool is_mat_hetero() const;
	VCFOneParentPhased *divide_by_positions(
								std::size_t first, std::size_t last) const;
	
	void impute();
	
private:
	Haplotype clip_haplotype(std::size_t sample_id, int i) const;
	Haplotype clip_ref_haplotype(std::size_t sample_id, int i) const;
	std::vector<Haplotype> collect_haplotypes_from_parents() const;
	std::vector<Haplotype> collect_haplotype_from_refs() const;
	std::vector<HaplotypePair> impute_cM(
								const std::vector<HaplotypePair>& prev_haps);
	
public:
	static VCFOneParentPhased *create(const STRVEC& samples, bool is_mat_phased,
						const VCFSmall *merged_vcf, const VCFSmall *orig_vcf,
						const Map& gmap, const VCFSmall *ref_vcf);
	static void impute_in_thread(void *config);
	static void impute_in_parallel(
				const std::vector<VCFOneParentPhased *>& vcfs, int num_threads);
	static STRVEC collect_samples(
						const std::vector<VCFOneParentPhased *>& vcfs);
	static VCFSmall *merge(const std::vector<VCFOneParentPhased *>& vcfs);
	static VCFSmall *impute_all(
				const std::vector<VCFOneParentPhased *>& vcfs, int num_threads);
};
#endif
