#include <iostream>

#include "../include/materials.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/common.h"

using namespace std;


//////////////////// Materials ////////////////////

Materials::Materials(const string& path_m, const Map *m,
						const string& path_p, const PedigreeTable *p) :
									path_map(path_m), geno_map(m),
									chr_maps(Map::create_chr_maps(m)),
									path_ped(path_p), ped(p) { }

Materials::~Materials() {
	delete geno_map;
	Common::delete_all(chr_maps);
	delete ped;
}

const Map *Materials::get_chr_map(int i) const {
	if(geno_map->is_empty())
		return chr_maps[0];
	else
		return chr_maps[i];
}

double Materials::total_cM() const {
	double	length = 0.0;
	for(auto p = chr_maps.begin(); p != chr_maps.end(); ++p)
		length += (*p)->total_cM();
	return length;
}

void Materials::display_map_info() const {
	cout << "Genetic Map : ";
	if(geno_map->is_empty()) {
		cout << "default map(1Mbp=1cM)." << endl;
	}
	else {
		cout << path_map << endl;
		cout << chr_maps.size() << " chrmosomes "
								<< total_cM() << " cM." << endl;
	}
}

Materials *Materials::create(const string& path_map, const string& path_ped) {
	const auto	*geno_map = Map::read(path_map);
	try {
		const auto	*ped = PedigreeTable::read(path_ped);
		return new Materials(path_map, geno_map, path_ped, ped);
	}
	catch(const std::exception& e) {
		delete geno_map;
		throw;
	}
}
