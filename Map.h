#ifndef __MAP
#define __MAP

#include <vector>


//////////////////// Map ////////////////////

class Map {
	class Record {
	public:
		std::string	chr;
		double		cM;
		double		Mbp;
		
	public:
		Record(const std::string& chr_, double cM_, double Mbp_) :
										chr(chr_), cM(cM_), Mbp(Mbp_) { }
		
		Record *copy() const { return new Record(chr, cM, Mbp); }
		
	public:
		static Record *create(const std::vector<std::string>& v);
	};
	
public:
	typedef std::vector<const Record *>::const_iterator	RIT;
	
private:
	const std::vector<const Record *>	records;
	
public:
	Map(const std::vector<const Record *>& rs) : records(rs) { }
	~Map();
	
	const std::vector<const Map *> divide_into_chromosomes() const;
	
	double bp_to_cM(long long bp) const;
	
private:
	RIT binary_search(RIT first, RIT last, double Mbp) const;
	
public:
	static Map *read(const std::string& path);
};
#endif
