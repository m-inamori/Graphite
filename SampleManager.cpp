#include <algorithm>
#include "SampleManager.h"
#include "common.h"

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
	this->imputed_samples.insert(samples.begin(), samples.end());
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
	// 片親だけimputed
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
		if(!this->is_parent_imputed_and_progenies_not_imputed(family))
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


vector<const Family *> SampleManager::make_families(const PedigreeTable *ped,
										const vector<string>& samples,
										const vector<size_t>& family_indices) {
	const vector<const Family *>	families = ped->make_families(samples);
	if(family_indices.empty())
		return families;
	
	vector<const Family *>	new_families;
	auto p = family_indices.begin();
	for(size_t i = 0; i < families.size(); ++i) {
		if(p == family_indices.end()) {
			;
		}
		else if(i != *p) {
			delete families[i];
		}
		else {
			new_families.push_back(families[i]);
			++p;
		}
	}
	return new_families;
}

SampleManager *SampleManager::create(const string& path_ped,
									const vector<string>& samples,
									int lower_progs,
									const vector<size_t>& family_indices) {
	const PedigreeTable	*ped_ = PedigreeTable::read(path_ped);
	const PedigreeTable	*ped = ped_->limit_samples(samples);
	delete ped_;
	
	const vector<const Family *>	families = make_families(ped, samples,
																family_indices);
	
	vector<const Family *>	large_families;
	vector<const Family *>	small_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family *family	= *p;
		if((int)family->num_progenies() >= lower_progs)
			large_families.push_back(family);
		else
			small_families.push_back(family);
	}
	
	return new SampleManager(ped, large_families, small_families);
}
