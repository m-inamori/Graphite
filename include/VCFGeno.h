// VCFGeno.h
#ifndef __VCFGENO
#define __VCFGENO

#include "GenoRecord.h"

class VCFGeno;


//////////////////// VCFGenoBase ////////////////////

class VCFGenoBase {
protected:
	const STRVEC	samples;
	const VCFSmall	*vcf;
	
public:
	VCFGenoBase(const STRVEC& s, const VCFSmall *v) : samples(s), vcf(v) { }
	virtual ~VCFGenoBase() { }
	
	virtual std::size_t size() const = 0;
	virtual GenoRecord *get_record(std::size_t i) const = 0;
	
	ll get_pos(std::size_t i) const { return get_record(i)->get_pos(); }
	const VCFSmall *get_ref_vcf() const { return vcf; }
	const STRVEC& get_samples() const { return samples; }
	std::size_t num_samples() const { return samples.size(); }
	
	std::vector<std::size_t> extract_columns(const STRVEC& samples) const;
	std::vector<int> extract_sample_genotypes(std::size_t sample_index) const;
	std::vector<int> clip_raw_haplotype(std::size_t sample_index,
														int side) const;
	void write(std::ostream& os, bool with_header=true) const;
	VCFGeno *extract_by_samples(const STRVEC& samples) const;
};


//////////////////// VCFGeno ////////////////////

class VCFGeno : public VCFGenoBase {
protected:
	std::vector<GenoRecord *>	records;
	
public:
	VCFGeno(const STRVEC& samples, const std::vector<GenoRecord *>& rs,
													const VCFSmall *vcf) :
									VCFGenoBase(samples, vcf), records(rs) { }
	~VCFGeno();
	
	///// virtual methods for VCFGenoBase /////
	GenoRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	std::size_t size() const override { return records.size(); }
	
	///// non-virtual methods /////
	const std::vector<GenoRecord *>& get_records() const { return records; }
	void clear_records() { records.clear(); }
	const VCFSmall *get_ref_vcf() const { return vcf; }
	
	static VCFGeno *convert(const VCFSmall *vcf);
	static VCFGeno *extract_samples(const STRVEC& samples, const VCFSmall *vcf);
	
	// join VCFs in order of samples
	static VCFGeno *join(const std::vector<const VCFGenoBase *>& vcfs,
													const STRVEC& samples);
	static VCFGeno *join(const VCFGenoBase *vcf1, const VCFGenoBase *vcf2,
														const STRVEC& samples);
	
	// Create a new VCF by obtaining the Genotypes of the parents from vcf1
	// and the progenies from vcf2
	static VCFGeno *create_by_two_vcfs(const VCFGenoBase *vcf1,
											const VCFSmall *vcf2,
											const STRVEC& samples);
};
#endif
