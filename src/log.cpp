#include <cassert>
#include <numeric>
#include <cmath>
#include "../include/log.h"

using namespace std;


//////////////////// Log ////////////////////

// Using log for simplicit
double Log::add(double log_a, double log_b) {
	// e.g. log_a = log2 log_b = log3
	// log(2 + 3) = log3 + log(2/3 + 1) = log3 + log(exp(log2-log3) + 1)
	if(log_a < log_b)
		return log_b + log(exp(log_a - log_b) + 1.0);
	else
		return log_a + log(exp(log_b - log_a) + 1.0);
}

double Log::sub(double log_a, double log_b) {
	// e.g. log_a = log3 log_b = log2
	// log(3 - 2) = log3 + log(1 - 2/3) = log3 + log(1 - exp(log2-log3))
	assert(log_a >= log_b);
	if(log_a == log_b || exp(log_b - log_a) == 1.0)
		return Log::LOGZERO;
	
	return log_a + log(1.0 - exp(log_b - log_a));
}

double Log::sum(const vector<double>& v) {
	if(v.empty())
		return Log::LOGZERO;
	return std::accumulate(v.begin() + 1, v.end(), v.front(), add);
}

double Log::diff(double log_a, double log_b) {
	if(log_a > log_b)
		return sub(log_a, log_b);
	else
		return sub(log_b, log_a);
}

map<string,double> Log::mat_log(const map<string,double>& M) {
	map<string,double>	result;
	for(auto p = M.begin(); p != M.end(); ++p) {
		const string&	k = p->first;
		double			v = p->second;
		result.insert(pair<string,double>(k, v > 0.0 ? log(v) : Log::LOGZERO));
	}
	return result;
}

double Log::modified_log(double x) {
	return x == 0.0 ? Log::LOGZERO : log(x);
}
