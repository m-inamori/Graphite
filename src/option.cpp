#include <iostream>
#include <stdexcept>
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// Option ////////////////////

bool Option::is_efficient_chrom(int i) const {
	if(chroms.empty())
		return true;
	else
		return std::find(chroms.begin(), chroms.end(), i) != chroms.end();
}

void Option::print_info() const {
	// required
	cerr << "input VCF : " << path_vcf << endl;
	cerr << "pedigree : " << path_ped << endl;
	cerr << "output VCF : " << path_out << endl;
	
	// optional
	cerr << "number of threads : " << num_threads << endl;
	cerr << "number of progenies for large family : " << lower_progs << endl;
	
	if(!imputes_isolated_samples) {
		cerr << "isolate samples will not be imputed." << endl;
		if(outputs_unimputed_samples)
			cerr << "but, outputs these samples" << endl;
	}
}

string Option::flag_value(const string& s, int argc, char **argv) {
	for(size_t i = 1; i < (size_t)(argc - 1); ++i) {
		if(argv[i] == s)
			return argv[i+1];
	}
	return string();
}

string Option::flags_value(const string& short_key, const string& long_key,
														int argc, char **argv) {
	const string	value = Option::flag_value(long_key, argc, argv);
	if(!value.empty()) {
		return value;
	}
	return Option::flag_value(short_key, argc, argv);
}

bool Option::exists(const string& s, int argc, char **argv) {
	for(size_t i = 3U; i < (size_t)argc; ++i) {
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
			throw std::runtime_error(f + " is invalid format.");
		}
	}
	return array;
}

vector<size_t> Option::get_families(int argc, char **argv) {
	const string	s = flags_value("-f", "--family", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

vector<size_t> Option::get_chroms(int argc, char **argv) {
	const string	s = flags_value("-c", "--chrom", argc, argv);
	if(s.empty())
		return vector<size_t>();
	else
		return parse_array(s);
}

int Option::get_num_threads(int argc, char **argv) {
	const string	s = flags_value("-t", "--num-threads", argc, argv);
	if(s.empty())
		return 1;
	else
		return stoi(s);
}

size_t Option::get_lower_progenies(int argc, char ** argv) {
	const string	s = flag_value("--lower-progs", argc, argv);
	if(s.empty())
		return 10;
	
	try {
		return stoul(s);
	}
	catch(const std::invalid_argument& e) {
		throw std::runtime_error(s + " is not an integer.");
	}
}

double Option::get_precision_ratio(int argc, char **argv) {
	const bool	b1 = exists("--fast", argc, argv);
	const bool	b2 = exists("--precision", argc, argv);
	const bool	b3 = exists("--precision-ratio", argc, argv);
	int	counter = 0;
	counter += b1 ? 1 : 0;
	counter += b2 ? 1 : 0;
	counter += b3 ? 1 : 0;
	if(counter >= 2) {
		throw std::runtime_error("--fast, --precision, and --precision cannot "
												"specified at the same time.");
	}
	else if(counter == 0) {
		return 1.0;
	}
	else if(b1) {
		return 0.1;
	}
	else if(b2) {
		return 10.0;
	}
	else {
		return stof(flag_value("--precision-ratio", argc, argv));
	}
}

Option *Option::create(int argc, char **argv) {
	try {
		// The top three are required arguments
		const string	path_vcf = flags_value("-i", "--input", argc, argv);
		if(path_vcf.empty())
			throw std::runtime_error("input VCF not specified.");
		
		const string	path_ped = flags_value("-p", "--pedigree", argc, argv);
		if(path_ped.empty())
			throw std::runtime_error("pedigree file not specified.");
		
		const string	path_out = flags_value("-o", "--output", argc, argv);
		if(path_out.empty())
			throw std::runtime_error("output VCF not specified.");
		
		// Optional
		const string	path_ref = flags_value("-r", "--ref", argc, argv);
		const string	path_map = flags_value("-m", "--map", argc, argv);
		const vector<size_t>	families = get_families(argc, argv);
		const vector<size_t>	chroms = get_chroms(argc, argv);
		const int	num_threads = get_num_threads(argc, argv);
		const size_t	lower_progs = get_lower_progenies(argc, argv);
		const double	prec_ratio = get_precision_ratio(argc,argv);
		const bool	impute_isolated = !exists("--not-impute-isolated",
																argc, argv);
		const bool	correct_isolated = !exists("--not-correct-isolated",
																argc, argv);
		const bool	out_isolated = exists("--out-isolated", argc, argv);
		if(impute_isolated && out_isolated)
			return NULL;
		
		return new Option(path_vcf, path_ref, path_ped, path_map, families,
							chroms, num_threads, lower_progs, prec_ratio,
							impute_isolated, correct_isolated,
							out_isolated, path_out);
	}
	catch(const std::runtime_error& e) {
		cerr << "error : " << e.what() << endl;
		return NULL;
	}
	catch(const std::invalid_argument& e) {
		return NULL;
	}
};

void Option::usage(char **argv) {
	cerr << "Usage:" << endl;
	cerr << "  graphite [options]" << endl;
	cerr << "" << endl;
	cerr << "  -i <path>/--input <path>     Input VCF file" << endl;
	cerr << "  -p <path>/--pedigree <path>  Pedigree file" << endl;
	cerr << "  -o <path>/--output <path>    Output file" << endl;
	cerr << "" << endl;
	cerr << "Options:" << endl;
	cerr << "  -m <path>/--map <path>       Input map file" << endl;
	cerr << "  -r <path>/--ref <path>       Reference VCF file" << endl;
	cerr << "  -t <int>/--num-threads <int> Number of threads (1)" << endl;
	cerr << "  --lower-progs <int>          Lower number of progenies" << endl;
	cerr << "                               considered to a large family" << endl;
	cerr << "  --not-impute-isolated        Samples that are isolated in the" << endl;
	cerr << "                               pedigree are not imputed" << endl;
	cerr << "  --not-correct-isolated       Not correct wrong genotypes of" << endl;
	cerr << "                               isolated samples" << endl;
	cerr << "  --out-isolated               Do not output isolated samples" << endl;
	cerr << "                               only valid when used in conjunction" << endl;
	cerr << "                               with --not-impute-isolated" << endl;
	cerr << "  -c <int>/--chrom <int>       Impute only the chromosome represented" << endl;
	cerr << "                               by the specified index (0-based)." << endl;
	cerr << "  --precision-ratio <float>    Control the runtime for small pedigree" << endl;
	cerr << "                               HMM analysis. Default is 1.0." << endl;
	cerr << "                               Larger values increase runtime." << endl;
	cerr << "  --fast                       Shortcut for setting --precision-ratio" << endl;
	cerr << "                               to 0.1. Optimized for faster runtime" << endl;
	cerr << "                               with reduced precision." << endl;
	cerr << "  --precision                  Shortcut for setting --precision-ratio" << endl;
	cerr << "                               to 10.0. Enhanced precision at the" << endl;
	cerr << "                               cost of increased runtime." << endl;
}
