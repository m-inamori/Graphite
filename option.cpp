#include <iostream>
#include "option.h"
#include "common.h"

using namespace std;


//////////////////// Option ////////////////////

bool Option::is_efficient_chrom(int i) const {
	if(chroms.empty())
		return true;
	else
		return std::find(chroms.begin(), chroms.end(), i) != chroms.end();
}

string Option::flag_value(const string& s, int argc, char **argv) {
	for(size_t i = 4U; i < (size_t)(argc - 1); ++i) {
		if(argv[i] == s)
			return argv[i+1];
	}
	return string();
}

bool Option::exists(const string& s, int argc, char **argv) {
	for(size_t i = 4U; i < (size_t)(argc - 1); ++i) {
		if(argv[i] == s)
			return true;
	}
	return false;
}

vector<size_t> Option::parse_array(const string& f) {
	vector<size_t>	array;
	const vector<string>	v = Common::split(f, ',');
	for(auto p = v.begin(); p != v.end(); ++p) {
		const vector<string>	w = Common::split(*p, ':');
		if(w.size() == 1U) {
			array.push_back(stoi(w[0]));
		}
		else if(w.size() == 2U) {
			for(size_t i = stoi(w[0]); i < (size_t)stoi(w[1]); ++i)
				array.push_back(i);
		}
		else {
			throw std::invalid_argument("");
		}
	}
	return array;
}

vector<size_t> Option::get_families(int argc, char **argv) {
	const string	s = flag_value("-f", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

vector<size_t> Option::get_chroms(int argc, char **argv) {
	const string	s = flag_value("-c", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

int Option::get_num_threads(int argc, char **argv) {
	const string	s = flag_value("-t", argc, argv);
	if(s.empty())
		return 1;
	else
		return stoi(s);
}

Option *Option::create(int argc, char **argv) {
	if(argc < 5 || 13 < argc)
		return NULL;
	
	try {
		const vector<size_t>	families = get_families(argc, argv);
		const vector<size_t>	chroms = get_chroms(argc, argv);
		const int	num_threads = get_num_threads(argc, argv);
		const bool	all_out = exists("-a", argc, argv);
		const bool	integ = exists("-i", argc, argv);
		return new Option(argv[1], argv[2], argv[3], families,
							chroms, num_threads, all_out, integ, argv[argc-1]);
	}
	catch(std::invalid_argument& e) {
		return NULL;
	}
};

void Option::usage(char **argv) {
	cerr << argv[0] << " VCF ped map [-t num_threads] "
				<< "[-f family indices] [-c chrom indices] "
				<< "[-a] [-i] out." << endl;
	cerr << "family indices: (index|first:last)[,(index|first:last),[..]]"
																	<< endl;
	cerr << "chrom indices: same as family indices." << endl;
	cerr << "-a: all records out." << endl;
	cerr << "-i: integrate same samples." << endl;
}
