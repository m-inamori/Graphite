#ifndef __VCFHMM
#define __VCFHMM

#include "Map.h"


//////////////////// VCFHMM ////////////////////

template<typename R>
class VCFHMM : public VCFMeasurable {
public:
	// (log of probability, previous hidden state)
	using DP = std::vector<std::pair<double, int>>;
	
protected:
	const std::vector<R *>&	records;
	const double	E[4][4];			// exhaust probabilities
	
public:
	VCFHMM(const std::vector<R *>& rs,
								const Map& map_, double w);
	~VCFHMM() { }
	
protected:
	std::vector<int> trace_back(const std::vector<DP>& dps) const;
	
	// genetic distance between two records
	double dist(const R *r1, const R *r2) const {
		const double	d = (cM(r2->get_pos()) - cM(r1->get_pos())) / 100;
		if(d != 0.0)
			return d;
		else	// probably outside map
			return (r2->get_pos() - r1->get_pos()) * 1e-6;
	}
};
#endif
