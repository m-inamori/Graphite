#include <cassert>
#include "common.h"
#include "Pedigree.h"
#include "VCF.h"

using namespace std;


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

vector<string> PedigreeTable::get_family_children(const string& mat,
												  const string& pat) const {
	vector<string>	children;
	pair<string,string>	parents(mat, pat);
	for(auto p = table.begin(); p != table.end(); ++p) {
		if((*p)->parents() == parents)
			children.push_back((*p)->get_name());
	}
	return children;
}

vector<pair<string,string>> PedigreeTable::extract_families() const {
	set<pair<string,string>>	s;
	for(auto p = table.begin(); p != table.end(); ++p)
		s.insert((*p)->parents());
	vector<pair<string,string>>	families(s.begin(), s.end());
	std::sort(families.begin(), families.end());
	return families;
}

const PedigreeTable *PedigreeTable::read(const string& path) {
	const auto	table = Common::read_csv(path, ' ');
	vector<const Progeny *>	progs;
	for(auto p = table.begin(); p != table.end(); ++p)
		progs.push_back(Progeny::create(*p));
	return new PedigreeTable(progs);
}

const PedigreeTable *PedigreeTable::create(const string& path, VCFBase *vcf) {
	const PedigreeTable	*ped = PedigreeTable::read(path);
	const PedigreeTable	*ped2 = ped->filter_with_parents(vcf->get_samples());
	
	map<pair<string,string>,int>	counter;
	for(auto p = ped2->table.begin(); p != ped2->table.end(); ++p) {
		const Progeny	*prog = *p;
		counter[prog->parents()] += 1;
	}
	
	set<pair<string,string>>	families;
	for(auto p = counter.begin(); p != counter.end(); ++p) {
		if(p->second >= 10)
			families.insert(p->first);
	}
	
	vector<const Progeny *>	progs;
	for(auto p = ped2->table.begin(); p != ped2->table.end(); ++p) {
		const Progeny	*prog = *p;
		if(families.find(prog->parents()) != families.end())
			progs.push_back(prog->copy());
	}
	
	delete ped2;
	delete ped;
	return new PedigreeTable(progs);
}
