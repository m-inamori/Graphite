#ifndef __OPTIONSMALL
#define __OPTIONSMALL

#include <vector>
#include <string>


//////////////////// OptionSmall ////////////////////

class OptionSmall {
public:
	const Map&		map;
	const int		num_threads;
	const double	precision_ratio;
	const bool		imputes_isolated_samples;
	const bool		outputs_unimputed_samples;
	
public:
	OptionSmall(const Map& m, int num_t, double pr, bool ii, bool ou) :
											map(m),
											num_threads(num_t),
											precision_ratio(pr),
											imputes_isolated_samples(ii),
											outputs_unimputed_samples(ou) { }
};
#endif
