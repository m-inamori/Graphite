#ifndef __VCFHMMREF
#define __VCFHMMREF

#include <array>

#include "GenoRecord.h"
#include "Map.h"


//////////////////// VCFHMMRef ////////////////////

class VCFHMMRef : public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	
protected:
	const std::array<std::array<double, 5>, 4>	E;	// exhaust probabilities
	
public:
	VCFHMMRef(const Map& map_, double w);
	~VCFHMMRef() { }
	
protected:
	std::vector<int> trace_back(const std::vector<DP>& dps) const;
	
	// genetic distance between two records
	double dist(const GenoRecord *r1, const GenoRecord *r2) const {
		const double	d = (cM(r2->get_pos()) - cM(r1->get_pos())) / 100;
		if(d != 0.0)
			return d;
		else	// probably outside map
			return (r2->get_pos() - r1->get_pos()) * 1e-6;
	}
};
#endif
