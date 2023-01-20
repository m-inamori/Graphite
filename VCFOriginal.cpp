#include <algorithm>
#include "VCFOriginal.h"
#include "Pedigree.h"

using namespace std;

STRVEC VCFOriginal::select_last_header_line(const VCFRecord *record) const {
	const STRVEC&	new_samples = record->get_samples();
	STRVEC	v(new_samples.size() + 9);
	std::copy(header.back().begin(), header.back().begin() + 9, v.begin());
	std::copy(new_samples.begin(), new_samples.end(), v.begin() + 9);
	return v;
}

vector<STRVEC> VCFOriginal::select_header(const VCFRecord *record) const {
	vector<STRVEC>	new_header = header;
	new_header.back() = select_last_header_line(record);
	return new_header;
}

vector<int> VCFOriginal::select_columns(const pair<string,string>& parents,
										const PedigreeTable& pedigree) const {
	const string&	mat = parents.first;
	const string&	pat = parents.second;
	const int	mat_column = this->find_column(mat);
	const int	pat_column = this->find_column(pat);
	const auto	progenies = pedigree.get_progenies(mat, pat);
	set<string>	set_progenies;
	for(auto p = progenies.begin(); p != progenies.end(); ++p)
		set_progenies.insert((*p)->get_name());
	
	vector<int>	columns = { mat_column, pat_column };
	int	i = 0;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		const string&	sample = *p;
		if(set_progenies.find(sample) != set_progenies.end())
			columns.push_back(i + 9);
		++i;
	}
	return columns;
}

vector<vector<int>> VCFOriginal::collect_family_columns(
								const vector<pair<string,string>>& families,
								const PedigreeTable& pedigree) const {
	vector<vector<int>>	family_columns;
	for(auto p = families.begin(); p != families.end(); ++p)
		family_columns.push_back(this->select_columns(*p, pedigree));
	return family_columns;
}

VCFOriginal *VCFOriginal::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>&	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	return new VCFOriginal(header, samples, reader);
}
