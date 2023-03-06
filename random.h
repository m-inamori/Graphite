// random.h

#ifndef INCLUDE_RANDOM
#define	INCLUDE_RANDOM

#include <vector>
#include <algorithm>

class Random {
	int	seed;
	
public:
	Random(int s=2) : seed(s) { }
	
	// https://en.wikipedia.org/wiki/Linear_congruential_generator
	void next() {
		seed = (int)(48271LL * seed % 0x7FFFFFFF);
	}
	
	double random() {
		next();
		return float(seed) / 0x7FFFFFFF;
	}
	
	int randint(int first, int last) {
		next();
		return first + seed % (last - first);
	}
	
	bool randbool() {
		next();
		return (randint(0, 1024) >> 9) == 1;
	}
	
	// choose between first and last(exclusive) without allowing duplication
	std::vector<int> sample(int first, int last, int n) {
		std::vector<int>	result;
		while((int)result.size() < n) {
			const int	e = randint(first, last);
			if(std::find(result.begin(), result.end(), e) == result.end())
				result.push_back(e);
		}
		return result;
	}
	
	template<typename T>
	std::vector<T> shuffle(const std::vector<T>& v) {
		// tree‚ðŽg‚¦‚Î‘¬‚­‚È‚é
		std::vector<T>	remains = v;
		std::vector<T>	shuffled;
		while(shuffled.size() < v.size()) {
			const int	i = randint(0, (int)remains.size());
			shuffled.push_back(remains[i]);
			remains.erase(remains.begin() + i);
		}
		return shuffled;
	}
};

#endif
