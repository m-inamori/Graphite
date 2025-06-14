#ifndef __SAMPLEMANAGER
#define __SAMPLEMANAGER

#include <set>
#include "KnownFamily.h"


//////////////////// SampleException ////////////////////

class SampleException : public ExceptionWithCode {
private:
    std::string	message;
	
public:
    SampleException() : message("error : large family not found") { }
    
    ErrorCode::Type get_error_code() const override {
		return ErrorCode::LARGE_FAMILY_NOT_FOUND;
	}
    const char *what() const noexcept override;
};


//////////////////// SampleManager ////////////////////

class SampleManager {
	const PedigreeTable	*ped;
	const std::vector<const KnownFamily *>	large_families;
	const std::vector<const KnownFamily *>	small_families;
	const std::size_t	lower_progs;
	std::set<std::string>	imputed_samples;
	
public:
	SampleManager(const PedigreeTable *p,
								const std::vector<const KnownFamily *>& f1,
								const std::vector<const KnownFamily *>& f2,
								std::size_t lower_p) :
									ped(p), large_families(f1),
									small_families(f2), lower_progs(lower_p) { }
	~SampleManager();
	
	const std::vector<const KnownFamily *>& get_large_families() const {
		return large_families;
	}
	const std::vector<const KnownFamily *>& get_small_families() const {
		return small_families;
	}
	
	bool is_imputed(const std::string& sample) const;
	bool is_all_progenies_imputed(const KnownFamily *family) const;
	bool is_known(const std::string& sample) const {
		return sample != "0";
	}
	bool is_unknown(const std::string& sample) const {
		return sample == "0";
	}
	bool is_parents_imputed_and_progenies_not_imputed(
												const Family *family) const;
	bool is_parent_imputed_and_progenies_not_imputed(
												const Family *family) const;
	bool is_progeny_imputed(const Family *family) const;
	std::vector<std::string> collect_reference() const;
	std::vector<const Progeny *> extract_unimputed_progenies(
												const Family * family) const;
	// families in which parents are phased and progenies are not phased
	std::vector<const KnownFamily *> extract_both_imputed_families() const;
	std::vector<const KnownFamily *> extract_imputed_and_known_families() const;
	std::vector<const KnownFamily *> extract_both_known_families() const;
	std::vector<const KnownFamily *> extract_one_imputed_families() const;
	std::vector<const KnownFamily *> extract_one_known_parent_families() const;
	std::vector<const KnownFamily *> extract_progenies_phased_families() const;
	bool is_all_not_imputed(const std::vector<std::string>& samples) const;
	std::vector<std::string> extract_isolated_samples() const;
	std::vector<std::string> extract_non_imputed_samples() const;
	std::vector<std::string> collect_large_family_parents() const;
	void display_info() const;
	
	void add_imputed_samples(const std::vector<std::string>& samples);
	void clear_imputed_samples();
	
public:
	static std::vector<const KnownFamily *> make_families(
								const PedigreeTable *ped,
								const std::vector<std::string>& samples,
								std::size_t lower_progs,
								const std::vector<std::size_t>& family_indices);
	static SampleManager *create(const PedigreeTable *ped,
								const std::vector<std::string>& samples,
								size_t lower_progs,
								const std::vector<std::size_t>& family_indices);
};
#endif
