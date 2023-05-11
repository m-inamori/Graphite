#ifndef __BAUMWELCHWITHFIXEDTS
#define __BAUMWELCHWITHFIXEDTS

#include <vector>
#include <string>
#include <map>


//////////////////// BaumWelch ////////////////////

class BaumWelch {
	typedef std::map<std::pair<char,char>,double>	Matrix;
	typedef std::vector<std::map<char,double>>	Table;
	
	const std::vector<char>&	hidden_states;
	const std::vector<char>&	states;
	const std::string&			seq;
	const std::vector<Matrix>	Ts;
	
	std::map<char,double>			pi;
//	Matrix	T;			// Transition Matrix
	Matrix	E;			// Emission Matirx
	
public:
	BaumWelch(const std::string& seq_,
				const std::vector<char>& hidden_states_,
				const std::vector<char>& states_,
				const std::vector<Matrix>& Ts_);
	
	Matrix reverse_probs(const Matrix& M) const;
	std::string Viterbi(const Matrix& A) const;
	
private:
	std::map<char,double> initialize_pi();
	Matrix initialize_transition_matrix();
	Matrix initialize_emission_matrix();
	bool is_valid_Emission_Matrix(const Matrix& E) const;
	Table compute_alpha(const Matrix& A) const;
	Table compute_beta(const Matrix& A) const;
	Table compute_gamma(const Table& a, const Table& b) const;
	Table compute_parameters(const Matrix& A);
	void update_probabilities(std::map<char,double>& pi1, Matrix& E1);
	bool is_near_parameters(const Matrix& E1,
							const std::map<char,double>& pi1,
							const Matrix& E2,
							const std::map<char,double>& pi2) const;
	bool converge_pi_E();
	bool converge_pi();
	std::map<char,double> update_pi(std::map<char,double>& pi);
	Matrix reverse_emission_matrix() const;

public:
	static std::string impute(const std::string& seq,
								const std::vector<char> hidden_states,
								const std::vector<char> states,
								const std::vector<Matrix>& Ts);
	
	static std::vector<Matrix> modify_log(const std::vector<Matrix>& Ms);
	
	static double get(const Matrix& M, char s1, char s2);
	static void set(Matrix& M, char s1, char s2, double v);
};
#endif
