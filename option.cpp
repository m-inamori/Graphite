#include <iostream>
#include "option.h"

using namespace std;


//////////////////// Option ////////////////////

Option *Option::create(int argc, char **argv) {
	if(argc < 5 || 8 < argc)
		return NULL;
	
	const bool	debug = string(argv[argc-2]) == "-b";
	int	num_threads;
	if(argc <= 6) {
		num_threads = 1;
		if(argc - (debug ? 1 : 0) != 5)
			return NULL;
	}
	else {
		if(string(argv[4]) != "-t")
			return NULL;
		if(argc - (debug ? 1 : 0) != 7)
			return NULL;
		try {
			num_threads = atoi(argv[5]);
		}
		catch(std::invalid_argument& e) {
			return NULL;
		}
	}
	
	return new Option(argv[1], argv[2], argv[3],
						num_threads, argv[argc-1], debug);
};

void Option::usage(char **argv) {
	cerr << argv[0] << " VCF ped map [-t num_threads] [-b] out." << endl;
}
