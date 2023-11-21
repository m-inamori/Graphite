#include <iostream>
#include <algorithm>
#include "../include/SampleManager.h"
#include "../include/common.h"

using namespace std;


//////////////////// SampleManager ////////////////////

SampleManager::~SampleManager() {
	delete ped;
	Common::delete_all(large_families);
	Common::delete_all(small_families);
}

bool SampleManager::is_imputed(const string& sample) const {
	return this->imputed_samples.find(sample) != this->imputed_samples.end();
}

const Family *SampleManager::get_large_family(
				const pair<std::string, std::string>& parents) const {
	for(auto p = large_families.begin(); p != large_families.end(); ++p) {
		const Family	*family = *p;
		if(family->parents() == parents)
			return family;
	}
	return NULL;
}

void SampleManager::add_imputed_samples(const vector<string>& samples) {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p != "0")
			this->imputed_samples.insert(*p);
	}
}

void SampleManager::clear_imputed_samples() {
	this->imputed_samples.clear();
}

bool SampleManager::is_parents_imputed_and_progenies_not_imputed(
												const Family *family) const {
	if(!this->is_imputed(family->get_mat()) ||
				!this->is_imputed(family->get_pat()))
		return false;
	
	const vector<const Progeny *>&	progenies = family->get_progenies();
	for(auto p = progenies.begin(); p != progenies.end(); ++p) {
		if(!this->is_imputed((*p)->get_name()))
			return true;
	}
	return false;
}

bool SampleManager::is_parent_imputed_and_progenies_not_imputed(
												const Family *family) const {
	// only one parent is imputed
	if(!(this->is_imputed(family->get_mat()) ^
				this->is_imputed(family->get_pat())))
		return false;
	
	const vector<const Progeny *>&	progenies = family->get_progenies();
	for(auto p = progenies.begin(); p != progenies.end(); ++p) {
		if(!this->is_imputed((*p)->get_name()))
			return true;
	}
	return false;
}

bool SampleManager::is_progeny_imputed(const Family *family) const {
	if(is_unknown(family->get_mat()) && is_unknown(family->get_pat()))
		return false;
	
	if(this->is_imputed(family->get_mat()) ||
					this->is_imputed(family->get_pat()))
		return false;
	
	const vector<const Progeny *>&	progenies = family->get_progenies();
	for(auto p = progenies.begin(); p != progenies.end(); ++p) {
		if(this->is_imputed((*p)->get_name()))
			return true;
	}
	return false;
}

vector<const Family *> SampleManager::extract_small_families() const {
	vector<const Family *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const Family	*family = *p;
		if(!this->is_parents_imputed_and_progenies_not_imputed(family))
			continue;
		
		vector<const Progeny *>	new_progenies;
		const vector<const Progeny *>&	progenies = family->get_progenies();
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			if(!this->is_imputed((*q)->get_name()))
				new_progenies.push_back((*q)->copy());
		}
		families.push_back(new Family(family->get_mat(),
										family->get_pat(), new_progenies));
	}
	return families;
}

vector<const Family *>
			SampleManager::extract_single_parent_phased_families() const {
	vector<const Family *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const Family	*family = *p;
		if(!this->is_parent_imputed_and_progenies_not_imputed(family) ||
									this->is_unknown(family->get_mat()) ||
									this->is_unknown(family->get_pat()))
			continue;
		
		vector<const Progeny *>	new_progenies;
		const vector<const Progeny *>&	progenies = family->get_progenies();
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			if(!this->is_imputed((*q)->get_name()))
				new_progenies.push_back((*q)->copy());
		}
		families.push_back(new Family(family->get_mat(),
										family->get_pat(), new_progenies));
	}
	return families;
}

// family in which one parent is phased and the other is unknown
vector<const Family *>
			SampleManager::extract_phased_and_unknown_parents_family() const {
	vector<const Family *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const Family	*family = *p;
		if(!this->is_parent_imputed_and_progenies_not_imputed(family) ||
									(this->is_known(family->get_mat()) &&
									 this->is_known(family->get_pat())))
			continue;
		
		vector<const Progeny *>	new_progenies;
		const vector<const Progeny *>&	progenies = family->get_progenies();
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			if(!this->is_imputed((*q)->get_name()))
				new_progenies.push_back((*q)->copy());
		}
		families.push_back(new Family(family->get_mat(),
										family->get_pat(), new_progenies));
	}
	return families;
}

vector<const Family *>
			SampleManager::extract_progenies_phased_families() const {
	vector<const Family *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const Family	*family = *p;
		if(!this->is_progeny_imputed(family))
			continue;
		
		vector<const Progeny *>	new_progenies;
		const vector<const Progeny *>&	progenies = family->get_progenies();
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			new_progenies.push_back((*q)->copy());
		}
		families.push_back(new Family(family->get_mat(),
										family->get_pat(), new_progenies));
	}
	return families;
}

bool SampleManager::is_all_not_imputed(const vector<string>& samples) const {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(is_imputed(*p))
			return false;
	}
	return true;
}

vector<string> SampleManager::extract_isolated_samples() const {
	// If there is a connected sample in the family
	// but the entire sample of the pedigree is not imputed,
	// it is considered isolated.
	vector<string>	samples;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const Family	*family = *p;
		const auto&	f_samples = family->get_samples();
		if(is_all_not_imputed(f_samples)) {
			for(auto q = f_samples.begin(); q != f_samples.end(); ++q) {
				if(is_known(*q) && !is_imputed(*q))
					samples.push_back(*q);
			}
		}
		else if(is_unknown(family->get_mat()) &&
							is_unknown(family->get_pat())) {
			const auto&	progs = family->get_progenies();
			for(auto q = progs.begin(); q != progs.end(); ++q) {
				if(!is_imputed((*q)->get_name()))
					samples.push_back((*q)->get_name());
			}
		}
	}
	return samples;
}

vector<string> SampleManager::get_large_parents() const {
	set<string>	s;
	for(auto p = large_families.begin(); p != large_families.end(); ++p) {
		const Family	*family = *p;
		s.insert(family->get_mat());
		s.insert(family->get_pat());
	}
	return vector<string>(s.begin(), s.end());
}

vector<string> SampleManager::collect_large_family_parents() const {
	set<string>	samples;
	for(auto p = large_families.begin(); p != large_families.end(); ++p) {
		samples.insert((*p)->get_mat());
		samples.insert((*p)->get_pat());
	}
	return vector<string>(samples.begin(), samples.end());
}

void SampleManager::display_info() const {
	cerr << ped->size() << " samples" << endl;
	
	if(large_families.size() == 1) {
		cerr << "1 large family";
	}
	else {
		cerr << large_families.size() << " large families";
	}
	cerr << " (number of progenies >= " << lower_progs << ")" << endl;
	
	if(small_families.size() == 1) {
		cerr << "1 small family" << endl;
	}
	else {
		cerr << small_families.size() << " small families" << endl;
	}
}

vector<const Family *> SampleManager::make_families(const PedigreeTable *ped,
										const vector<string>& samples,
										const vector<size_t>& family_indices) {
	const vector<const Family *>	families = ped->make_families(samples);
	
	vector<const Family *>	new_families;
	const set<string>	set_samples(samples.begin(), samples.end());
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		// Families without any progeny in the VCF are deleted.
		const auto&	progs = family->get_progenies();
		vector<const Progeny *>	filtered_progs;
		for(auto q = progs.begin(); q != progs.end(); ++q) {
			const Progeny	*prog = *q;
			if(set_samples.find(prog->get_name()) != set_samples.end())
				filtered_progs.push_back(prog->copy());
		}
		if(!filtered_progs.empty()) {
			// Treat parents who are not in the VCF as unknown
			const string	mat = set_samples.find(family->get_mat()) !=
									set_samples.end() ? family->get_mat() : "0";
			const string	pat = set_samples.find(family->get_pat()) !=
									set_samples.end() ? family->get_pat() : "0";
			new_families.push_back(new Family(mat, pat, filtered_progs));
		}
	}
	
	Common::delete_all(families);
	
	if(family_indices.empty())
		return new_families;
	
	vector<const Family *>	new_families2;
	auto p = family_indices.begin();
	for(size_t i = 0; i < new_families.size(); ++i) {
		if(p == family_indices.end()) {
			;
		}
		else if(i != *p) {
			delete new_families[i];
		}
		else {
			new_families2.push_back(new_families[i]);
			++p;
		}
	}
	return new_families2;
}

SampleManager *SampleManager::create(const string& path_ped,
									const vector<string>& samples,
									size_t lower_progs,
									const vector<size_t>& family_indices) {
	const PedigreeTable	*ped_ = PedigreeTable::read(path_ped);
	const PedigreeTable	*ped = ped_->limit_samples(samples);
	delete ped_;
	
	const vector<const Family *>	families = make_families(ped, samples,
																family_indices);
	
	const set<string>	set_samples(samples.begin(), samples.end());
	vector<const Family *>	large_families;
	vector<const Family *>	small_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family *family	= *p;
		if(family->num_progenies() >= lower_progs &&
					set_samples.find(family->get_mat()) != set_samples.end() &&
					set_samples.find(family->get_pat()) != set_samples.end())
			large_families.push_back(family);
		else
			small_families.push_back(family);
	}
	
	return new SampleManager(ped, large_families, small_families, lower_progs);
}
