#ifndef __OPTION
#define __OPTION

#include <vector>
#include <string>


//////////////////// Option ////////////////////

class Option {
public:
	const std::string	path_vcf;
	const std::string	path_ped;
	const std::string	path_map;
	const std::vector<std::size_t>	families;
	const int			chroms;
	const int			num_threads;
	const std::string	path_out;
	
public:
	Option(const std::string& vcf, const std::string& ped,
			const std::string& map, const std::vector<std::size_t>& fam,
			int chr, int nt, const std::string& out) :
				path_vcf(vcf), path_ped(ped), path_map(map),
				families(fam), chroms(chr), num_threads(nt), path_out(out) { }
public:
	static Option *create(int argc, char **argv);
	static void usage(char **argv);
	
private:
	static std::string flag_value(const std::string& s, int argc, char **argv);
	static std::vector<std::size_t> parse_families(const std::string& f);
	static std::vector<std::size_t> get_families(int argc, char **argv);
	static int get_chroms(int argc, char **argv);
	static int get_num_threads(int argc, char **argv);
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
