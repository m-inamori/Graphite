#ifndef __MAP
#define __MAP

#include <vector>

#include "exception_with_code.h"


//////////////////// Map ////////////////////

class Map {
public:
	class Record {
	public:
		std::string	chr;
		double		cM;
		double		Mbp;
		
	public:
		Record(const std::string& chr_, double cM_, double Mbp_) :
										chr(chr_), cM(cM_), Mbp(Mbp_) { }
		
		Record *copy() const { return new Record(chr, cM, Mbp); }
		
		std::string to_string() const;
		
	public:
		static Record *create(const std::vector<std::string>& v);
	};
	
public:
	typedef std::vector<const Record *>::const_iterator	RIT;
	
private:
	const std::vector<const Record *>	records;
	
public:
	explicit Map(const std::vector<const Record *>& rs) : records(rs) { }
	~Map();
	
	bool is_empty() const { return records.empty(); }
	double total_cM() const { return records.back()->cM; }
	const std::vector<const Map *> divide_into_chromosomes() const;
	
	void check_in_order() const;
	double bp_to_cM(long long bp) const;
	
private:
	RIT binary_search(RIT first, RIT last, double Mbp) const;
	
public:
	static std::vector<std::vector<std::string>> read_lines(
													const std::string& path);
	static Map *read(const std::string& path);
	static Map *default_map();
	static std::vector<const Map *> create_chr_maps(const Map *m);
	static double Kosambi(double d);
};


//////////////////// VCFMeasurable ////////////////////

class VCFMeasurable {
	const Map&	gmap;
	
public:
	explicit VCFMeasurable(const Map& m) : gmap(m) { }
	~VCFMeasurable() { }
	
	const Map& get_map() const { return gmap; }
	double cM(long long bp) const { return gmap.bp_to_cM(bp); }
	double total_cM() const { return gmap.total_cM(); }
};


//////////////////// MapFormatException ////////////////////

class MapFormatException : public ExceptionWithCode {
private:
	std::string	message;
	
public:
	explicit MapFormatException(const std::vector<std::string>& lines);
	
	ErrorCode::Type get_error_code() const override {
		return ErrorCode::MAP_INVALID_FORMAT;
	}
	
	const char *what() const noexcept override {
		return message.c_str();
	}
};


//////////////////// TwiceChrException ////////////////////

class TwiceChrException : public ExceptionWithCode {
private:
	std::string	message;
	
public:
	explicit TwiceChrException(const std::string& chr);
	
	ErrorCode::Type get_error_code() const override {
		return ErrorCode::MAP_INVALID_FORMAT;
	}
	
	const char *what() const noexcept override {
		return message.c_str();
	}
};


//////////////////// OutOfOrderCMException ////////////////////

class OutOfOrderCMException : public ExceptionWithCode {
private:
	std::string	message;
	
public:
	OutOfOrderCMException(const Map::Record *r1, const Map::Record *r2);
	
	ErrorCode::Type get_error_code() const override {
		return ErrorCode::MAP_INVALID_FORMAT;
	}
	
	const char *what() const noexcept override {
		return message.c_str();
	}
};


//////////////////// OutOfOrderMbpException ////////////////////

class OutOfOrderMbpException : public ExceptionWithCode {
private:
	std::string	message;
	
public:
	OutOfOrderMbpException(const Map::Record *r1, const Map::Record *r2);
	
	ErrorCode::Type get_error_code() const override {
		return ErrorCode::MAP_INVALID_FORMAT;
	}
	
	const char *what() const noexcept override {
		return message.c_str();
	}
};

#endif
