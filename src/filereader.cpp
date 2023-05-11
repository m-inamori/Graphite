#include "../include/filereader.h"

using namespace std;


//////////////////// FileReaderBase ////////////////////

FileReaderBase *FileReaderBase::create(const string& path) {
	if(path.substr(path.length() - 3U) != ".gz")
		return new FileReader(path);
	else
		return new FileReaderGZ(path);
}


//////////////////// FileReader ////////////////////

FileReader::FileReader(const string& path_) : FileReaderBase(),
											path(path_), ifs(path_.c_str()) {
	if(!ifs) {
		cerr << "error : can't open " << path << "." << endl;
	}
	
	iter = buffer.begin();
}

bool FileReader::getline(std::string& line) {
	std::string	l;
	if(iter == buffer.end()) {
		buffer.clear();
		while(std::getline(ifs, l)) {
			buffer.push_back(l);
			if(buffer.size() == upper_lines)
				break;
		}
		iter = buffer.begin();
	}
	
	if(iter ==  buffer.end())
		return false;
	
	line = *iter;
	++iter;
	return true;
}


//////////////////// FileReaderGZ ////////////////////

FileReaderGZ::FileReaderGZ(const string& path_) :
								FileReaderBase(), path(path_) {
	fz = gzopen(path_.c_str(), "r");
	if(fz == NULL) {
		cerr << "error : can't open " << path << "." << endl;
	}
	
	iter = buffer.begin();
}

FileReaderGZ::~FileReaderGZ() {
	gzclose(fz);
}

bool FileReaderGZ::getline(string& line) {
	if(iter == buffer.end()) {
		buffer.clear();
		string	buff;
		while(getline_core(buff)) {
			if(buff.substr(buff.length()-1, 1).c_str()[0] == '\n')
				buff.pop_back();
			buffer.push_back(buff);
			if(buffer.size() == upper_lines)
				break;
		}
		iter = buffer.begin();
	}
	
	if(iter ==  buffer.end())
		return false;
	
	line = *iter;
	++iter;
	return true;
}

bool FileReaderGZ::getline_core(string& line) {
	char	buff[4096];
	line = "";
	while(true) {
		if(!gzgets(fz, buff, 4096))
			return !line.empty();
		line = line + buff;
		if(line.substr(line.length()-1, 1).c_str()[0] == '\n')
			return true;
	}
	return false;
}
