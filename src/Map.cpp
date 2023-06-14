#include <stdexcept>
#include "../include/common.h"
#include "../include/Map.h"

using namespace std;


//////////////////// Map::Record ////////////////////

Map::Record *Map::Record::create(const vector<string>& v) {
	if(v.size() != 3)
		return NULL;
	
	try {
		const string&	chr_ = v[0];
		const double	cM_ = stod(v[1]);
		const double	Mbp_ = stod(v[2]);
		return new Record(chr_, cM_, Mbp_);
	}
	catch(std::invalid_argument& e) {
		return NULL;
	}
}


//////////////////// Map ////////////////////

Map::~Map() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

vector<const Map *> Map::create_chr_maps(const Map *m) {
	if(m->is_empty())
		// Make one map and use it for every chromosome.
		// That makes it easier to delete.
		return vector<const Map *>(1, Map::default_map());
	else
		return m->divide_into_chromosomes();
}

double Map::bp_to_cM(long long bp) const {
	const double	Mbp = bp * 1e-6;
	RIT	p = binary_search(records.begin(), records.end()-1, Mbp);
	const Record	*r1 = *p;
	const Record	*r2 = *(p + 1);
	if(r1->Mbp == r2->Mbp)
		return r1->cM;
	
	return r1->cM + (r2->cM - r1->cM) / (r2->Mbp - r1->Mbp) * (Mbp - r1->Mbp);
}

const vector<const Map *> Map::divide_into_chromosomes() const {
	vector<const Map *>	maps;
	string	prev_chr = "";
	vector<const Record *>	rs;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const string&	chr = (*p)->chr;
		if(chr != prev_chr) {
			if(!rs.empty()) {
				maps.push_back(new Map(rs));
				rs.clear();
			}
			prev_chr = chr;
		}
		rs.push_back((*p)->copy());
	}
	maps.push_back(new Map(rs));
	return maps;
}

Map::RIT Map::binary_search(RIT first, RIT last, double Mbp) const {
	if(first == last - 1)
		return first;
	
	const RIT	mid = first + (last - first) / 2;
	if(Mbp < (*mid)->Mbp)
		return binary_search(first, mid, Mbp);
	else
		return binary_search(mid, last, Mbp);
}

Map *Map::read(const string& path) {
	vector<const Record *>	records;
	if(path.empty())
		return new Map(records);
	
	const auto	table = Common::read_csv(path);
	for(auto p = table.begin(); p != table.end(); ++p) {
		records.push_back(Record::create(*p));
	}
	return new Map(records);
}

Map *Map::default_map() {
	const Record	*r1 = new Record("1", 0.0, 0.0);
	const Record	*r2 = new Record("1", 1.0, 1.0);
	const vector<const Record *>	rs = { r1, r2 };
	return new Map(rs);
}
