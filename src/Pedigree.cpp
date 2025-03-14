#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cassert>
#include "../include/common.h"
#include "../include/Pedigree.h"

using namespace std;


//////////////////// FormatException ////////////////////

FormatException::FormatException(const vector<string>& lines) {
	stringstream	ss;
	if(lines.size() == 1)
		ss << "error : the following line doesn't have four columns :";
	else
		ss << "error : the following lines don't have four columns :";
	
	for(auto p = lines.begin(); p != lines.end(); ++p)
		ss << '\n' << *p;
	
	message = ss.str();
}

const char *FormatException::what() const noexcept {
	return message.c_str();
}


//////////////////// ParentsException ////////////////////

ParentsException::ParentsException(const vector<string>& parents) {
	stringstream	ss;
	if(parents.size() == 1)
		ss << "error : the following parent isn't defined :";
	else
		ss << "error : the following parents aren't defined :";
	
	for(auto p = parents.begin(); p != parents.end(); ++p)
		ss << '\n' << *p;
	
	message = ss.str();
}

const char *ParentsException::what() const noexcept {
	return message.c_str();
}


//////////////////// SamplesException ////////////////////

SamplesException::SamplesException(const vector<string>& samples) {
	stringstream	ss;
	if(samples.size() == 1)
		ss << "error : the following sample isn't in pedigree :";
	else
		ss << "error : the following samples aren't in pedigree :";
	
	for(auto p = samples.begin(); p != samples.end(); ++p)
		ss << '\n' << *p;
	
	message = ss.str();
}

const char *SamplesException::what() const noexcept {
	return message.c_str();
}


//////////////////// Family ////////////////////

Family::~Family() {
	for(auto p = progenies.begin(); p != progenies.end(); ++p)
		delete *p;
}

vector<string> Family::collect_samples() const {
	vector<string>	samples;
	samples.push_back(this->mat);
	samples.push_back(this->pat);
	for(auto p = this->progenies.begin(); p != this->progenies.end(); ++p)
		samples.push_back((*p)->get_name());
	return samples;
}


//////////////////// Progeny ////////////////////

const Progeny *Progeny::copy() const {
	return new Progeny(this->family, this->name, this->mat, this->pat);
}

bool Progeny::with_parents() const {
	return mat != "0" && pat != "0";
}

bool Progeny::is_all_in_samples(const set<string>& samples) const {
	return samples.find(name) != samples.end() &&
		   samples.find(mat) != samples.end() &&
		   samples.find(pat) != samples.end();
}

const Progeny *Progeny::create(const vector<string>& v) {
	assert(v.size() >= 4U);
	return new Progeny(v[0], v[1], v[2], v[3]);
}


//////////////////// PedigreeTable ////////////////////

PedigreeTable::~PedigreeTable() {
	for(auto p = table.begin(); p != table.end(); ++p)
		delete *p;
}

const PedigreeTable *PedigreeTable::filter_with_parents(
									const vector<string>& samples) const {
	vector<const Progeny *>	progs;
	set<string>	set_samples(samples.begin(), samples.end());
	for(auto p = table.begin(); p != table.end(); ++p) {
		if((*p)->with_parents() && (*p)->is_all_in_samples(set_samples))
			progs.push_back((*p)->copy());
	}
	// not destruct Progeny in destructor
	return new PedigreeTable(progs);
}

vector<const Progeny *> PedigreeTable::get_progenies(const string& mat,
													  const string& pat) const {
	vector<const Progeny *>	children;
	set<string>	set_children;	// for avoiding duplication
	pair<string,string>	parents(mat, pat);
	for(auto p = table.begin(); p != table.end(); ++p) {
		if((*p)->parents() == parents) {
			if(set_children.find((*p)->get_name()) == set_children.end()) {
				children.push_back((*p)->copy());
				set_children.insert((*p)->get_name());
			}
		}
	}
	return children;
}

vector<const Family *> PedigreeTable::extract_families() const {
	// collect parents
	set<pair<string,string>>	s;
	for(auto p = this->table.begin(); p != this->table.end(); ++p)
		s.insert((*p)->parents());
	
	vector<const Family *>	families;
	for(auto p = s.begin(); p != s.end(); ++p) {
		const string&	mat = p->first;
		const string&	pat = p->second;
		const auto	progenies = this->get_progenies(mat, pat);
		const Family	*family = new Family(mat, pat, progenies);
		families.push_back(family);
	}
	return families;
}

vector<const Family *> PedigreeTable::make_families(
										const vector<string>& samples) const {
	map<pair<string, string>, vector<const Progeny *>>	progs;
	set<string>	used;
	for(auto p = this->table.begin(); p != this->table.end(); ++p) {
		const Progeny	*prog = *p;
		if(used.find(prog->get_name()) == used.end()) {
			progs[prog->parents()].push_back(prog->copy());
			used.insert(prog->get_name());
		}
	}
	
	vector<const Family *>	families;
	const set<string>	set_samples(samples.begin(), samples.end());
	for(auto p = progs.begin(); p != progs.end(); ++p) {
		const auto	parents = p->first;
		const auto	progenies = p->second;
		// Parent not found in samples of VCF is treated as unknown
		const string	mat = parents.first;
		const string	pat = parents.second;
		if(mat != "0" && pat != "0") {
			const Family	*family = new Family(mat, pat, progenies);
			families.push_back(family);
		}
		else {
			for(auto q = progenies.begin(); q != progenies.end(); ++q) {
				vector<const Progeny *>	new_progs{ *q };
				const Family	*family = new Family(mat, pat, new_progs);
				families.push_back(family);
			}
		}
	}
	std::sort(families.begin(), families.end(),
				[](const Family *lh, const Family *rh)
				{ return lh->parents() < rh->parents(); });
	return families;
}

vector<string> PedigreeTable::check_samples_in_pedigree(
								const vector<string>& samples) const {
	set<string>	ped_samples;
	for(auto p = table.begin(); p != table.end(); ++p)
		ped_samples.insert((*p)->get_name());
	vector<string>	missing_samples;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(ped_samples.find(*p) == ped_samples.end())
			missing_samples.push_back(*p);
	}
	return missing_samples;
}

const PedigreeTable *PedigreeTable::limit_samples(
								const vector<string>& samples) const {
	const auto	missing_samples = check_samples_in_pedigree(samples);
	if(!missing_samples.empty()) {
		throw SamplesException(missing_samples);
	}
	
	set<string>	set_samples(samples.begin(), samples.end());
	vector<const Progeny *>	progs;
	for(auto p = table.begin(); p != table.end(); ++p) {
		if(set_samples.find((*p)->get_name()) != set_samples.end())
			progs.push_back((*p)->copy());
	}
	return new PedigreeTable(progs);
}

vector<string> PedigreeTable::check_parents() const {
	set<string>	set_progs;
	for(auto p = table.begin(); p != table.end(); ++p)
		set_progs.insert((*p)->get_name());
	
	set<string>	missing_parents;
	for(auto p = table.begin(); p != table.end(); ++p) {
		const string&	mat = (*p)->get_mat();
		if(mat != "0" && set_progs.find(mat) == set_progs.end())
			missing_parents.insert(mat);
		const string&	pat = (*p)->get_pat();
		if(pat != "0" && set_progs.find(pat) == set_progs.end())
			missing_parents.insert(pat);
	}
	return vector<string>(missing_parents.begin(), missing_parents.end());
}

vector<vector<string>> PedigreeTable::read_lines(const string& path) {
	ifstream	ifs(path.c_str());
	if(!ifs)
		throw FileNotFoundException(path);
	
	vector<vector<string>>	table;
	vector<string>	not_four_columns_lines;
	string	line;
	while(getline(ifs, line)) {
		if(line.c_str()[line.length()-1] == '\r')
			line = line.substr(0, line.length()-1);
		istringstream	iss(line);
		vector<string>	v;
		string	token;
		while(iss >> token) {
			v.push_back(token);
		}
		if(v.size() != 4) {
			not_four_columns_lines.push_back(line);
		}
		table.push_back(v);
	}
	if(!not_four_columns_lines.empty()) {
		throw FormatException(not_four_columns_lines);
	}
	return table;
}

const vector<const Progeny *> PedigreeTable::remove_duplidated_progenies(
										const vector<const Progeny *>& progs) {
	vector<const Progeny *>	new_progs;
	set<string>	used_names;
	for(auto p = progs.begin(); p != progs.end(); ++p) {
		const string&	prog_name = (*p)->get_name();
		if(used_names.find(prog_name) == used_names.end()) {
			new_progs.push_back(*p);
			used_names.insert(prog_name);
		}
		else {
			cerr << "Progeny " << prog_name << " is duplicated." << endl;
			delete *p;
		}
	}
	return new_progs;
}

const PedigreeTable *PedigreeTable::create(
								const vector<const Progeny *>& progs) {
	const auto	ped = new PedigreeTable(progs);
	const auto	missing_parents = ped->check_parents();
	if(!missing_parents.empty()) {
		delete ped;
		throw ParentsException(missing_parents);
	}
	return ped;
}

const PedigreeTable *PedigreeTable::read(const string& path) {
	const auto	table = PedigreeTable::read_lines(path);
	vector<const Progeny *>	progs;
	for(auto p = table.begin(); p != table.end(); ++p)
		progs.push_back(Progeny::create(*p));
	const auto	unique_progs = remove_duplidated_progenies(progs);
	return create(unique_progs);
}
