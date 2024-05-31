#include <sstream>

#include "../include/exception_with_code.h"

using namespace std;


//////////////////// FileNotFoundException ////////////////////

FileNotFoundException::FileNotFoundException(const string& path) {
	stringstream	ss;
	
	ss << "error : " << path << " can't open.";
	message = ss.str();
}

const char *FileNotFoundException::what() const noexcept {
	return message.c_str();
}
