#ifndef __MATERIALS
#define __MATERIALS

#include <vector>
#include <string>
#include <exception>

class Map;
class PedigreeTable;


//////////////////// Materials ////////////////////

class Materials {
	const std::string	path_map;
	const Map	*geno_map;
	const std::vector<const Map *>	chr_maps;
	
	const std::string	path_ped;
	const PedigreeTable	*ped;
	
public:
	Materials(const std::string& path_m, const Map *m,
				const std::string& path_p, const PedigreeTable *p);
	~Materials();
	
	const std::string& get_path_map() const { return path_map; }
	const Map& get_map() const { return *geno_map; }
	const Map *get_chr_map(int i) const;
	double total_cM() const;
	void display_map_info() const;
	const std::string& get_path_ped() const { return path_ped; }
	const PedigreeTable *get_ped() const { return ped; }
	
public:
	static Materials *create(const std::string& path_map,
								const std::string& path_ped);
};

#endif
