#include <cmath>
#include "../include/VCFHMM.h"
#include "../include/VCFFamily.h"

using namespace std;

template<typename R>
VCFHMM<R>::VCFHMM(const std::vector<R *>& rs,
									const Map& map_, double w) :
				VCFMeasurable(map_), records(rs),
				E{{log(1.0-w*2), log(w/2),     log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(1.0-w*2), log(w/2),     log(w)},
				  {log(w/2),     log(w/2),     log(1.0-w*2), log(w)}} { }

template<typename R>
vector<int> VCFHMM<R>::trace_back(const vector<DP>& dps) const {
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
	hs[M-2] = prev_h;
	for(int i = (int)M - 2; i > 0; --i) {
		prev_h = dps[i][prev_h].second;
		hs[i-1] = prev_h;
	}
	return hs;
}

template class VCFHMM<VCFRecord>;
template class VCFHMM<VCFFamilyRecord>;
