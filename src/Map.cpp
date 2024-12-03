#include <fstream>
#include <sstream>
#include <cmath>
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
	catch(const std::invalid_argument& e) {
		vector<string>	lines(1, Common::join(v, ','));
		throw MapFormatException(lines);
	}
}

string Map::Record::to_string() const {
	stringstream	ss;
	ss << chr << ',' << cM << ',' << Mbp;
	return ss.str();
}


//////////////////// Map ////////////////////

Map::~Map() {
	for(auto p = records.begin(); p != records.end(); ++p) {
		if(*p != NULL)
			delete *p;
	}
}

vector<const Map *> Map::create_chr_maps(const Map *m) {
	if(m->is_empty())
		// Make one map and use it for every chromosome.
		// That makes it easier to delete.
		return vector<const Map *>(1, Map::default_map());
	else
		return m->divide_into_chromosomes();
}

void Map::check_in_order() const {
	string	prev_chr = records[0]->chr;
	double	prev_cM = records[0]->cM;
	double	prev_Mbp = records[0]->Mbp;
	vector<string>	past_chrs;
	for(size_t i = 1; i < records.size(); ++i) {
		const string	chr = records[i]->chr;
		const double	cM = records[i]->cM;
		const double	Mbp = records[i]->Mbp;
		if(prev_chr != chr) {
			past_chrs.push_back(prev_chr);
			if(std::find(past_chrs.begin(), past_chrs.end(), chr)
													!= past_chrs.end()) {
				throw TwiceChrException(chr);
			}
			prev_chr = chr;
		}
		else {
			if(cM < prev_cM)
				throw OutOfOrderCMException(records[i-1], records[i]);
			else if(Mbp < prev_Mbp)
				throw OutOfOrderMbpException(records[i-1], records[i]);
		}
		prev_cM = cM;
		prev_Mbp = Mbp;
	}
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

vector<vector<string>> Map::read_lines(const string& path) {
	ifstream	ifs(path.c_str());
	if(!ifs)
		throw FileNotFoundException(path);
	
	vector<vector<string>>	table;
	vector<string>	not_three_columns_lines;
	string	line;
	while(getline(ifs, line)) {
		if(line.c_str()[line.length()-1] == '\r')
			line = line.substr(0, line.length()-1);
		istringstream	iss(line);
		vector<string>	v = Common::split(line, ',');
		if(v.size() != 3 && not_three_columns_lines.size() < 5) {
			not_three_columns_lines.push_back(line);
		}
		table.push_back(v);
	}
	if(!not_three_columns_lines.empty()) {
		throw MapFormatException(not_three_columns_lines);
	}
	return table;
}

Map *Map::read(const string& path) {
	vector<const Record *>	records;
	if(path.empty())
		return new Map(records);
	
	const auto	table = read_lines(path);
	try {
		for(auto p = table.begin(); p != table.end(); ++p) {
			records.push_back(Record::create(*p));
		}
		if(records.size() < 2) {
			vector<string>	lines;
			for(auto p = records.begin(); p != records.end(); ++p)
				lines.push_back((*p)->to_string());
			throw MapFormatException(lines);
		}
	}
	catch(const ExceptionWithCode& e) {
		Common::delete_all(records);
		throw;
	}
	auto	*geno_map = new Map(records);
	try {
		geno_map->check_in_order();
	}
	catch(const ExceptionWithCode& e) {
		delete geno_map;
		throw;
	}
	return geno_map;
}

Map *Map::default_map() {
	const Record	*r1 = new Record("1", 0.0, 0.0);
	const Record	*r2 = new Record("1", 1.0, 1.0);
	const vector<const Record *>	rs = { r1, r2 };
	return new Map(rs);
}

double Map::Kosambi(double d) {
	return tanh(2*d) / 2;
}


//////////////////// MapFormatException ////////////////////

MapFormatException::MapFormatException(const vector<string>& lines) {
	stringstream	ss;
	if(lines.size() == 1)
		ss << "error : the following line doesn't have three columns :";
	else
		ss << "error : the following lines don't have three columns :";
	
	for(auto p = lines.begin(); p != lines.end(); ++p)
		ss << '\n' << *p;
	
	message = ss.str();
}


//////////////////// TwiceChrException ////////////////////

TwiceChrException::TwiceChrException(const string& chr) {
	stringstream	ss;
	ss << "error : the following chromosome comes out twice\n";
	ss << chr;
	message = ss.str();
}


//////////////////// OutOfOrderCMException ////////////////////

OutOfOrderCMException::OutOfOrderCMException(const Map::Record *r1,
											 const Map::Record *r2) {
	stringstream	ss;
	ss << "error : the following cMs are out of order\n";
	ss << r1->to_string() << '\n';
	ss << r2->to_string();
	message = ss.str();
}


//////////////////// OutOfOrderMbpException ////////////////////

OutOfOrderMbpException::OutOfOrderMbpException(const Map::Record *r1,
											   const Map::Record *r2) {
	stringstream	ss;
	ss << "error : the following Mbps are out of order\n";
	ss << r1->to_string() << '\n';
	ss << r2->to_string();
	message = ss.str();
}
