#include <cmath>

#include "../include/VCFHMMRef.h"

using namespace std;


//////////////////// VCFHMMRef ////////////////////

VCFHMMRef::VCFHMMRef(const Map& map_, double w) :
	VCFMeasurable(map_),
	E{{
	{{log(1.0-w*2), log(w/3),         log(w/3),         log(w/3),     log(w)}},
	{{log(w/3),     log(0.9-5.3/3*w), log(0.1+0.1*w),   log(w/3),     log(w)}},
	{{log(w/3),     log(0.1+0.1*w),   log(0.9-5.3/3*w), log(w/3),     log(w)}},
	{{log(w/3),     log(w/3),         log(w/3),         log(1.0-w*2), log(w)}}
    }}
	{ }

vector<int> VCFHMMRef::trace_back(const vector<VCFHMMRef::DP>& dps) const {
	const size_t	M = dps.size();
	vector<int>	hs(M, 0);
	
	pair<double, int>	max_pair(-1e300, 0);
	int	max_h = 0;
	const size_t	L = dps[0].size();
	for(int h = 0; h < (int)L; ++h) {
		if(dps.back()[h] > max_pair) {
			max_pair = dps.back()[h];
			max_h = h;
		}
	}
	hs[M-1] = max_h;
	int	prev_h = dps.back()[max_h].second;
	if(M >= 2) {
		hs[M-2] = prev_h;
		for(int i = (int)M - 2; i > 0; --i) {
			prev_h = dps[i][prev_h].second;
			hs[i-1] = prev_h;
		}
	}
	return hs;
}
