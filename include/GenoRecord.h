// GenoRecord.h
#ifndef __GENORECORD
#define __GENORECORD

#include "VCF.h"
#include "Genotype.h"


//////////////////// GenoRecord ////////////////////

class GenoRecord {
protected:
	const ll	pos;
	std::vector<int>	geno;
	
public:
	GenoRecord(ll p, const std::vector<int>& g) : pos(p), geno(g) { }
	virtual ~GenoRecord() { }
	
	const ll get_pos() const { return pos; }
	const std::vector<int>& get_geno() const { return geno; }
	std::size_t num_samples() const { return geno.size(); }
	int unphased(std::size_t i) const { return Genotype::unphased(geno[i]); }
	bool is_hetero(std::size_t i) const { return Genotype::is_hetero(geno[i]); }
	bool is_ref_homo(std::size_t i) const {
		return Genotype::is_ref_homo(geno[i]);
	}
	bool is_alt_homo(std::size_t i) const {
		return Genotype::is_alt_homo(geno[i]);
	}
	bool is_homo(std::size_t i) const { return Genotype::is_homo(geno[i]); }
	bool is_NA(std::size_t i) const { return Genotype::is_NA(geno[i]); }
	bool is_00(std::size_t i) const { return Genotype::is_00(geno[i]); }
	bool is_01(std::size_t i) const { return Genotype::is_01(geno[i]); }
	bool is_11(std::size_t i) const { return Genotype::is_11(geno[i]); }
	bool is_phased(std::size_t i) const { return Genotype::is_phased(geno[i]); }
	
	// allele of j(0 or 1) side
	int get_allele(std::size_t i, int j) const {
		return Genotype::get_allele(geno[i], j);
	}
	
	std::vector<int> unphased_gts() const;
	
	void set_geno(std::size_t i, int gt) { geno[i] = gt; }
	
	void write(const VCFRecord *record, std::ostream& os) const;
};
#endif
