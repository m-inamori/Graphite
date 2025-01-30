#ifndef __VCFHMM
#define __VCFHMM

#include "Map.h"

class VCFRecord;
class VCFFamilyRecord;


//////////////////// VCFHMM ////////////////////

class VCFHMM : public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	
protected:
	const std::vector<VCFFamilyRecord *>&	records;
	const double	E[4][4];			// exhaust probabilities
	
public:
	VCFHMM(const std::vector<VCFFamilyRecord *>& rs,
								const Map& map_, double w);
	~VCFHMM() { }
	
protected:
	std::vector<int> trace_back(const std::vector<DP>& dps) const;
	
	// genetic distance between two records
	double dist(const VCFRecord *r1, const VCFRecord *r2) const;
};
#endif
