#ifndef __VCFORIGINAL
#define __VCFORIGINAL

#include "VCF.h"

class PedigreeTable;


//////////////////// VCFOriginal ////////////////////

class VCFOriginal : public VCFHuge {
public:
	VCFOriginal(const std::vector<STRVEC>& h, const STRVEC& s, VCFReader *r) :
															VCFHuge(h, s, r) { }
	
	std::vector<STRVEC> select_header(const VCFRecord *record) const;
	std::vector<int> select_columns(
							const std::pair<std::string,std::string>& parents,
							const PedigreeTable& pedigree) const;
	std::vector<std::vector<int>> collect_family_columns(
				const std::vector<std::pair<std::string,std::string>>& families,
				const PedigreeTable& pedigree) const;
	
private:
	STRVEC select_last_header_line(const VCFRecord *record) const;
	
public:
	static VCFOriginal *read(const std::string& path);
};
#endif
