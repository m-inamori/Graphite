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
	const std::vector<std::size_t>	chroms;
	const int			num_threads;
	const std::string	path_out;
	const double		ratio;
	const std::size_t	lower_progs;
	const bool			only_large_families;
	
public:
	Option(const std::string& vcf, const std::string& ped,
			const std::string& map, const std::vector<std::size_t>& fam,
			const std::vector<std::size_t>& chr,
			int nt, std::size_t lp, bool ol, const std::string& out) :
									path_vcf(vcf), path_ped(ped), path_map(map),
									families(fam), chroms(chr), num_threads(nt),
									path_out(out), ratio(0.01), lower_progs(lp),
									only_large_families(ol) { }
	
	bool is_efficient_chrom(int i) const;
	
public:
	static Option *create(int argc, char **argv);
	static void usage(char **argv);
	
private:
	static std::string flag_value(const std::string& s, int argc, char **argv);
	static bool exists(const std::string& s, int argc, char **argv);
	static std::vector<std::size_t> parse_array(const std::string& f);
	static std::vector<std::size_t> get_families(int argc, char **argv);
	static std::vector<std::size_t> get_chroms(int argc, char **argv);
	static int get_num_threads(int argc, char **argv);
	static std::size_t get_lower_progenies(int argc, char ** argv);
};


//////////////////// OptionImpute ////////////////////

class OptionImpute {
public:
	const int	max_dist;
	const int	min_positions;
	const int	min_graph;
	const double	min_crossover;
	
public:
	OptionImpute(int md, int mp, int mg, double mc) :
				max_dist(md), min_positions(mp),
				min_graph(mg), min_crossover(mc) { }
};
#endif
