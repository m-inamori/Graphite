#ifndef __MATHCOMMON
#define __MATHCOMMON

#include <functional>


//////////////////// MathCommon ////////////////////

namespace MathCommon {
	template<class F>
	double Simpson(F f, double a, double b, int n) {
		double s = f(a) + f(b);
		
		for (int k = 1; k < 2 * n; k += 2) {
			s += 4 * f((k * a + (2 * n - k) * b) / (2 * n));
		}
		
		for (int k = 2; k < 2 * n; k += 2) {
			s += 2 * f((k * a + (2 * n - k) * b) / (2 * n));
		}
		
	    return s * (b - a) / (2 * n) / 3;
	}
}
#endif
