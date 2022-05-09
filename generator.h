#ifndef __GENERATOR
#define __GENERATOR

// emulator of Python's generator
template <typename T>
class Generator {
public:
	Generator() { }
	virtual ~Generator() { }
	
	virtual int size() const = 0;
	virtual const T *next() = 0;
};
#endif
