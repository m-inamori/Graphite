#ifndef __BIASPROBABILITY
#define __BIASPROBABILITY

#include <vector>
#include <map>


//////////////////// BiasProbability ////////////////////

class BiasProbability {
	typedef std::map<std::pair<int,int>,double>	LineProbs;
	
	const double	probability;
	std::map<std::pair<int,int>,double>	memo_C;		// log(nCk)
	// { num samples: [{ (bias, max bias): probability }] }
	std::map<int,std::vector<LineProbs>>	ps;
	
public:
	BiasProbability(double p) : probability(p) { }
	~BiasProbability() { }
	
	int compute_max_bias(int N, double float_cM);
	
private:
	// log(nCk)
	double log_C(int n, int k);
	LineProbs init(int N);
	void extend(std::vector<LineProbs>& subps, int cM, int N);
	int compute_max_bias_core(const LineProbs& lps);
	LineProbs step(const LineProbs& lps, int N);
	
public:
	static int bias(int i, int N);
	
private:
	static BiasProbability	*bias_probability;
};
#endif
