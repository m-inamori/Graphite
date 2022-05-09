#ifndef __FILEREADER
#define __FILEREADER

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <zlib.h>
#include "generator.h"


//////////////////// FileReaderBase ////////////////////

class FileReaderBase {
public:
	FileReaderBase() { }
	virtual ~FileReaderBase() { }
	
	virtual bool is_opened() const = 0;
	virtual bool getline(std::string& line) = 0;
	virtual const std::string& get_path() const = 0;
	
public:
	static FileReaderBase *create(const std::string& path);
};


//////////////////// FileReader ////////////////////

class FileReader : public FileReaderBase {
	static const size_t	upper_lines = 1000U;
	const std::string	path;	// for debug
	std::ifstream	ifs;
	std::vector<std::string>	buffer;
	std::vector<std::string>::const_iterator	iter;
	
public:
	FileReader(const std::string& path_);
	~FileReader() { ifs.close(); }
	
	bool is_opened() const { return (bool)ifs; }
	bool getline(std::string& line);
	const std::string& get_path() const { return path; }
};


//////////////////// FileReaderGZ ////////////////////

class FileReaderGZ : public FileReaderBase {
	static const size_t	upper_lines = 1000U;
	const std::string	path;	// for debug
	gzFile	fz;
	std::vector<std::string>	buffer;
	std::vector<std::string>::const_iterator	iter;
	
public:
	FileReaderGZ(const std::string& path_);
	~FileReaderGZ();
	
	bool is_opened() const { return fz != NULL; }
	bool getline(std::string& line);
	const std::string& get_path() const { return path; }
	
private:
	bool getline_core(std::string& line);
};
#endif
