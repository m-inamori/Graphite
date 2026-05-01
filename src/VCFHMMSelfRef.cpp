#include <cmath>

#include "../include/VCFHMMSelfRef.h"
#include "../include/VCFFamily.h"

using namespace std;


//////////////////// VCFHMMSelfRef ////////////////////

VCFHMMSelfRef::VCFHMMSelfRef(const Map& map_, double w) :
										VCFMeasurable(map_),
										E(calc_E(w))
										{ }

array<array<double, 5>, 4> VCFHMMSelfRef::calc_E(double w) {
	// hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
	// observed 0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3, N/A: 4
	// Assume a 0.1 probability of confusing 0|1 and 1|0.
	array<array<double, 5>, 4>	E1 {};
	for(int gt = 0; gt < 4; ++gt) {
		array<double, 5>	orig_prob;
		for(int k = 0; k < 4; ++k) {
			const int	i = k >> 1;
			const int	j = k & 1;
			const int	a1 = (gt >> i) & 1;
			const int	a2 = (gt >> j) & 1;
			const int	prog_gt = a1 | (a2 << 1);
			if(prog_gt == 1) {
				orig_prob[1] += 0.25 * 0.9;
				orig_prob[2] += 0.25 * 0.1;
			}
			else if(prog_gt == 2) {
				orig_prob[1] += 0.25 * 0.1;
				orig_prob[2] += 0.25 * 0.9;
			}
			else {
				orig_prob[gt] += 0.25;
			}
		}
		
		// put errors
		for(int prog_gt = 0; prog_gt < 4; ++prog_gt) {
			for(int error_gt = 0; error_gt < 5; ++error_gt) {
				if(error_gt == 4)
					E1[gt][error_gt] = w;
				else if(error_gt == prog_gt)
					E1[gt][error_gt] += orig_prob[gt] * (1.0 - w*2);
				else
					E1[gt][error_gt] += orig_prob[gt] * (w/3);
			}
		}
	}
	
	array<array<double, 5>, 4> E;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 5; ++j) {
				E[i][j] = log(E1[i][j]);
		}
	}
	return E;
}

vector<int> VCFHMMSelfRef::trace_back(const vector<DP>& dps) const {
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
