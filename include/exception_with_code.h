#ifndef __EXCEPTIONWITHCODE
#define __EXCEPTIONWITHCODE

#include <exception>
#include <string>

#include "error_codes.h"


//////////////////// ExceptionWithCode ////////////////////

// with error code
class ExceptionWithCode : public std::exception {
public:
	ExceptionWithCode() { }
	
	virtual ErrorCode::Type get_error_code() const = 0;
	virtual const char* what() const noexcept { return ""; }
};


//////////////////// FileNotFoundException ////////////////////

class FileNotFoundException : public ExceptionWithCode {
private:
    std::string	message;
	
public:
    FileNotFoundException(const std::string& path);
    
    ErrorCode::Type get_error_code() const {
		return ErrorCode::FILE_NOT_FOUND;
	}
    const char *what() const noexcept override;
};

#endif
