#ifndef __VCFONEPARENTPHASED
#define __VCFONEPARENTPHASED

#include "VCFFamily.h"
#include "option.h"
#include "graph.h"

class Map;
class PedigreeTable;
class VCFOriginal;


//////////////////// VCFOneParentPhased ////////////////////

class VCFOneParentPhased : public VCFFamily {
	const bool	is_mat_phased;
	
protected:
	const Map&	genetic_map;
	
public:
	VCFOneParentPhased(const std::vector<STRVEC>& h, const STRVEC& s,
										std::vector<VCFFamilyRecord *> rs,
										bool is_mat_phased, const Map& m);
	~VCFOneParentPhased() { }
	
	bool is_mat_hetero() const;
	
	void impute();
	
private:
	double cM(std::size_t i) const;
	char determine_which_comes_from(VCFFamilyRecord *record,
												std::size_t i) const;
	std::string make_seq(std::size_t i) const;
	std::string impute_sample_seq(std::size_t i,
							const std::vector<double>& cMs, double min_c) const;
	void update_each(std::size_t i, std::size_t j, char c);
	void update(std::size_t i, const std::string& seq);
	
public:
	static VCFFamily *impute_by_parent(const VCFSmall *orig_vcf,
										const VCFSmall *parent_imputed_vcf,
										const STRVEC& samples,
										bool is_mat_phased, const Map& gmap);
	// samplesに対応するRecordを作る
	static VCFRecord *merge_records(std::vector<VCFFamily *>& vcfs,
										std::size_t i, const STRVEC& samples);
};
#endif
