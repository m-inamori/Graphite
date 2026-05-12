#ifndef __VCFIMPHETEROHOMO
#define __VCFIMPHETEROHOMO

#include "VCFImputable.h"
#include "VCFFillableRecord.h"
#include "Map.h"


//////////////////// VCFImpHeteroHomo ////////////////////

class VCFImpHeteroHomo : public VCFImputable, public VCFMeasurable {
	std::vector<VCFFillableRecord *>	records;
	bool	is_mat_hetero;
	
public:
	VCFImpHeteroHomo(const STRVEC& s,
						const std::vector<VCFFillableRecord *>& rs,
						bool is_mat_hetero_, const Map& m,
						const VCFSmall *vcf) :
									VCFImputable(s, vcf),
									VCFMeasurable(m),
									records(rs),
									is_mat_hetero(is_mat_hetero_) { }
	VCFImpHeteroHomo(const VCFImpHeteroHomo&) = delete;
	VCFImpHeteroHomo& operator=(const VCFImpHeteroHomo&) = delete;
	~VCFImpHeteroHomo() { }
	
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
	std::size_t imputed_index() const { return is_mat_hetero ? 0 : 1; }
	std::size_t non_imputed_index() const { return is_mat_hetero ? 1 : 0; }
	
	int update_each(std::size_t i, std::size_t j, char c) const;
	void update(std::size_t i, const std::vector<std::string>& seqs);
	char determine_haplotype(int which_zero,
								int homo_int_gt, int prog_int_gt) const;
	std::string make_seq(std::size_t i) const;
	std::string impute_sample_seq(std::size_t j, const std::vector<double>& cMs,
															double min_c) const;
	
	///// virtual methods for VCFImputable /////
	std::size_t amount() const override { return 1; }
	void impute() override;
};
#endif
