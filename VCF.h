#ifndef __VCF
#define __VCF

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "filereader.h"


typedef long long					ll;
typedef std::vector<std::string>	STRVEC;
typedef std::pair<int, ll>			POSITION;


//////////////////// VCFRecord ////////////////////

class VCFRecord {
protected:
	STRVEC	v;
	const STRVEC	samples;
	
public:
	VCFRecord(const STRVEC& v_, const STRVEC& s_) : v(v_), samples(s_) { }
	~VCFRecord();
	
	const STRVEC& get_v() const { return v; }
	const STRVEC& get_samples() const { return samples; }
	const std::string& chrom() const { return this->v[0]; }
	std::pair<std::string,ll>	position() const;
	ll	pos() const { return stoll(this->v[1]); }
	const std::string	format() const { return this->v[8]; }
	std::string	gt(const std::string& sample) const;
	STRVEC gts() const;
	std::string get_gt(std::size_t i) const { return v[i+9]; }
	std::string get_GT(std::size_t i) const { return v[i+9].substr(0, 3); }
	int get_int_gt(std::size_t i) const;
	std::vector<int> get_int_gts() const;
	STRVEC extract_v(const STRVEC& samples) const;
	void write(std::ostream& os) const;
	
	void copy_properties(STRVEC::iterator it) const;
	void set_GT(std::size_t i, const std::string& gt);
	void set_GTs(const STRVEC& GTs);
	void set(const STRVEC& new_v) { v = new_v; }
};


//////////////////// VCFBase ////////////////////

class VCFBase {
protected:
	const std::vector<STRVEC>	header;
	const STRVEC	samples;
	const std::map<std::string,std::size_t>	sample_ids;
	mutable std::map<std::string,int>	chrs;
	
public:
	// Pythonよりheaderからsamplesを取り出すのがめんどうなので
	VCFBase(const std::vector<STRVEC>& header_, const STRVEC& samples_);
	~VCFBase() { }
	
	const std::vector<STRVEC>& get_header() const { return header; }
	const STRVEC& get_samples() const { return samples; }
	std::map<std::string,std::size_t>
	number_samples(const STRVEC& samples_) const;
	POSITION position(const std::pair<std::string,ll>& p) const;
	POSITION record_position(const VCFRecord& record) const;
	std::string chr(int chr_id) const;
	int find_column(const std::string& sample) const;
	void write_header(std::ostream& os) const;
	
	void copy_chrs(VCFBase *vcf) { vcf->chrs = chrs; }
	
private:
	void determine_chromosome_id();
};


//////////////////// VCFReader ////////////////////

class VCFReader {
	FileReaderBase	*reader;
	std::vector<STRVEC>	header;
	
public:
	VCFReader(const std::string& path);
	~VCFReader();
	
	void read_header();
	const std::vector<STRVEC>& get_header() { return header; }
	STRVEC get_samples() const;
	STRVEC next();
};


//////////////////// VCFHuge ////////////////////

class VCFHuge : public VCFBase {
	VCFReader	*reader;
	
public:
	VCFHuge(const std::vector<STRVEC>& header_, const STRVEC& samples_,
													VCFReader *reader_);
	~VCFHuge();
	
	VCFRecord *next();
	VCFRecord *proceed(POSITION pos);
	
public:
	static VCFHuge *read(const std::string& path);
};


//////////////////// VCFSmall ////////////////////

class VCFSmall : public VCFBase {
protected:
	std::vector<VCFRecord *>	records;
	
public:
	VCFSmall(const std::vector<STRVEC>& header, const STRVEC& samples,
											std::vector<VCFRecord *> rs);
	~VCFSmall();
	
	const std::vector<VCFRecord *>& get_records() const { return records; }
	std::size_t size() const { return records.size(); }
	void write(std::ostream& os) const;
	
	void update_genotypes(const std::vector<STRVEC>& GT_table);
	
public:
	static VCFSmall *read(const std::string& path);
};
#endif
