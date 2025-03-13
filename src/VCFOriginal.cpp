#include <algorithm>
#include "../include/VCFOriginal.h"
#include "../include/Pedigree.h"

using namespace std;

VCFOriginal *VCFOriginal::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>&	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	return new VCFOriginal(header, samples, reader);
}
