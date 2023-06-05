#ifndef __VCFPROGENYPHASED
#define __VCFPROGENYPHASED

#include "VCFFamily.h"

class Family;


//////////////////// VCFOneParentPhased ////////////////////

class VCFProgenyPhased : public VCFFamily {
public:
	typedef std::vector<std::size_t>	PPI;
	
	struct ConfigThread {
		const VCFSmall	*orig_vcf;
		const VCFSmall	*merged_vcf;
		const std::vector<std::pair<const Family *, PPI>>&	families;
		const std::size_t	first;
		const int	num_threads;
		std::vector<VCFProgenyPhased *>&	results;
		
		ConfigThread(const VCFSmall *o, const VCFSmall *m,
					 const std::vector<std::pair<const Family *, PPI>>& fs,
					 std::size_t f, int T,
					 std::vector<VCFProgenyPhased *>& rs) :
			 						orig_vcf(o), merged_vcf(m),
			 						families(fs), first(f),
			 						num_threads(T), results(rs) { }
		
		std::size_t size() const { return families.size(); }
	};
	
private:
	const std::vector<std::size_t>	phased_progeny_indices;
	
public:
	VCFProgenyPhased(const std::vector<STRVEC>& h, const STRVEC& s,
										std::vector<VCFFamilyRecord *> rs,
										const std::vector<std::size_t>& ppi);
	~VCFProgenyPhased() { }
	
	void determine_parent(bool is_mat);
	void impute();
	
private:
	static void impute_in_thread(void *config);
	
public:
	static VCFProgenyPhased *impute_by_progeny(const VCFSmall *orig_vcf,
										const VCFSmall *imputed_vcf,
										const STRVEC& samples,
										const std::vector<std::size_t>& ppi);
	static std::vector<VCFProgenyPhased *> impute_all_by_progeny(
					const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
					const std::vector<std::pair<const Family *, PPI>>& families,
					int num_threads);
};
#endif
