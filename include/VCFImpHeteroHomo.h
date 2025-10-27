#ifndef __VCFIMPHETEROHOMO
#define __VCFIMPHETEROHOMO

#include "VCFHeteroHomoOnePhased.h"

class Map;


//////////////////// VCFImpHeteroHomo ////////////////////

class VCFImpHeteroHomo : public VCFHeteroHomoOnePhased {
public:
	VCFImpHeteroHomo(const STRVEC& s,
						const std::vector<VCFFillableRecord *>& rs,
						bool is_mat_hetero, const Map& m,
						const VCFSmall *vcf) :
				VCFHeteroHomoOnePhased(s, rs, is_mat_hetero, m, vcf) { }
	~VCFImpHeteroHomo() { }
	
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
	
	///// virtual methods /////
	void impute() override;
};
#endif
