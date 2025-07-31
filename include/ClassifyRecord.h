#ifndef __CLASSIFYRECORD
#define __CLASSIFYRECORD


#include <map>

#include "TypeDeterminer.h"


enum class WrongType {
	RIGHT, MODIFIABLE, UNMODIFIABLE, MIX, UNSPECIFIED
};

class VCFRecord;
class VCFFamilyRecord;



//////////////////// ClassifyRecord ////////////////////

class ClassifyRecord {
public:
	// probability and genotypes of parents
	typedef std::pair<double, ParentComb>	GTComb;
	
private:
	// Singleton
	static ClassifyRecord	*instance;
	
	std::map<std::pair<std::size_t, double>, const TypeDeterminer *>	memo;
	
private:
	ClassifyRecord() { }
	~ClassifyRecord();
	
public:
	static ClassifyRecord *get_instance();
	
	const TypeDeterminer *get_TypeDeterminer(std::size_t n, double alpha);
	std::pair<ParentComb, WrongType> classify(const VCFFamilyRecord *record,
													const TypeDeterminer *td,
													bool one_parent);
	std::pair<ParentComb, WrongType> classify_self_record(
												const VCFRecord *record,
												const TypeDeterminer *td);
	
private:
	std::array<int, 3> count_int_gts(const std::vector<int>& gts) const;
	std::array<int, 3> count_int_gts(std::vector<int>::const_iterator first,
								std::vector<int>::const_iterator last) const;
	std::vector<GTComb> filter_pairs(const std::vector<GTComb>& combs) const;
	WrongType select_wrong_type(ParentComb comb, int mat_gt,
								int pat_gt, bool one_parent) const;
	bool is_matched(int mat_gt, int pat_gt, ParentComb comb) const;
	std::pair<ParentComb, WrongType> select_pair(std::vector<GTComb>& combs,
												int mat_gt, int pat_gt) const;
	std::pair<ParentComb, WrongType> classify_record_core(
								const std::vector<GTComb>& combs_,
								int mat_gt, int pat_gt, bool one_parent) const;
};
#endif
