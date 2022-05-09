#ifndef __LOG
#define __LOG

#include <string>
#include <vector>
#include <map>


namespace Log {

//////////////////// Constants ////////////////////

const double	LOGZERO = -10000.0;


//////////////////// Log ////////////////////

double add(double log_a, double log_b);
double sub(double log_a, double log_b);
double sum(const std::vector<double>& v);
double diff(double log_a, double log_b);
std::map<std::string,double> mat_log(const std::map<std::string,double>& M);
double modified_log(double x);

}
#endif
