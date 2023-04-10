#ifndef __SAMPLEMANAGER
#define __SAMPLEMANAGER

#include <set>
#include "pedigree.h"


//////////////////// SampleManager ////////////////////

class SampleManager {
	const PedigreeTable	*ped;
	const std::vector<const Family *>	large_families;
	const std::vector<const Family *>	small_families;
	std::set<std::string>	imputed_samples;
	
public:
	SampleManager(const PedigreeTable *p, const std::vector<const Family *>& f1,
										const std::vector<const Family *>& f2):
							ped(p), large_families(f1), small_families(f2) { }
	~SampleManager();
	
	const std::vector<const Family *>& get_large_families() const {
		return large_families;
	}
	const std::vector<const Family *>& get_small_families() const {
		return small_families;
	}
	const Family *get_large_family(
				const std::pair<std::string, std::string>& parents) const;
	
	bool is_imputed(const std::string& sample) const;
	bool is_parents_imputed_and_progenies_not_imputed(
												const Family *family) const;
	bool is_parent_imputed_and_progenies_not_imputed(
												const Family *family) const;
	// 補完されていないが両親は補完されている家系
	std::vector<const Family *> extract_small_families() const;
	// 補完されていないが片方の親だけ補完されている家系
	std::vector<const Family *> extract_single_parent_phased_families() const;
	
	void add_imputed_samples(const std::vector<std::string>& samples);
	void clear_imputed_samples();
	
public:
	static std::vector<const Family *> make_families(const PedigreeTable *ped,
								const std::vector<std::string>& samples,
								const std::vector<std::size_t>& family_indices);
	static SampleManager *create(const std::string& path_ped,
								const std::vector<std::string>& samples,
								int lower_progs,
								const std::vector<std::size_t>& family_indices);
};
#endif