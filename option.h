#ifndef __OPTION
#define __OPTION

#include <string>


//////////////////// Option ////////////////////

class Option {
public:
	const std::string	path_vcf;
	const std::string	path_ped;
	const std::string	path_map;
	const int			num_threads;
	const std::string	path_out;
	const bool			debug;
	
public:
	Option(const std::string& vcf, const std::string& ped,
			const std::string& map, int nt,
			const std::string& out, bool d) :
				path_vcf(vcf), path_ped(ped), path_map(map),
				num_threads(nt), path_out(out), debug(d) { }
	public:
		static Option *create(int argc, char **argv);
		static void usage(char **argv);
};


//////////////////// OptionImpute ////////////////////

class OptionImpute {
public:
	const int	max_dist;
	const int	min_positions;
	const int	min_graph;
	
public:
	OptionImpute(int md, int mp, int mg) :
				max_dist(md), min_positions(mp), min_graph(mg) { }
};
#endif
