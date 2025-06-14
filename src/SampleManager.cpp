#include <iostream>
#include <algorithm>
#include "../include/SampleManager.h"
#include "../include/KnownFamily.h"
#include "../include/common.h"

using namespace std;


//////////////////// SampleException ////////////////////

const char *SampleException::what() const noexcept {
	return message.c_str();
}


//////////////////// SampleManager ////////////////////

SampleManager::~SampleManager() {
	Common::delete_all(large_families);
	Common::delete_all(small_families);
}

bool SampleManager::is_imputed(const string& sample) const {
	return this->imputed_samples.find(sample) != this->imputed_samples.end();
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

vector<string> SampleManager::collect_reference() const {
	set<string>	ref_samples;
	for(auto p = large_families.begin(); p != large_families.end(); ++p) {
		const auto	parents = (*p)->known_parents();
		ref_samples.insert(parents.begin(), parents.end());
	}
	return vector<string>(ref_samples.begin(), ref_samples.end());
}

vector<const Progeny *> SampleManager::extract_unimputed_progenies(
												const Family *family) const {
	vector<const Progeny *>	new_progenies;
	const vector<const Progeny *>&	progenies = family->get_progenies();
	for(auto q = progenies.begin(); q != progenies.end(); ++q) {
		if(!this->is_imputed((*q)->get_name()))
			new_progenies.push_back((*q)->copy());
	}
	return new_progenies;
}

vector<const KnownFamily *>
				SampleManager::extract_both_imputed_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		if(!this->is_parents_imputed_and_progenies_not_imputed(family))
			continue;
		
		families.push_back(family->create(extract_unimputed_progenies(family)));
	}
	return families;
}

vector<const KnownFamily *>
			SampleManager::extract_imputed_and_known_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		if(!this->is_parent_imputed_and_progenies_not_imputed(family) ||
									this->is_unknown(family->get_mat()) ||
									this->is_unknown(family->get_pat()))
			continue;
		
		families.push_back(family->create(extract_unimputed_progenies(family)));
	}
	return families;
}

vector<const KnownFamily *> SampleManager::extract_both_known_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		const auto	progs = extract_unimputed_progenies(family);
		if(this->is_imputed(family->get_mat()) ||
				this->is_imputed(family->get_pat()) ||
				this->is_unknown(family->get_mat()) ||
				this->is_unknown(family->get_pat()) ||
				progs.empty()) {
			Common::delete_all(progs);
			continue;
		}
		
		families.push_back(family->create(progs));
	}
	return families;
}

vector<const KnownFamily *>
				SampleManager::extract_one_imputed_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		if((this->is_imputed(family->get_mat()) &&
									!family->is_pat_known()) ||
				(this->is_imputed(family->get_pat()) &&
									!family->is_mat_known())) {
			families.push_back(family);
		}
	}
	
	vector<const KnownFamily *>	new_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*f = *p;
		const auto&	progs = extract_unimputed_progenies(f);
		for(auto q = progs.begin(); q != progs.end(); ++q) {
			const auto	*g = new KnownFamily(f->get_mat(), f->get_pat(),
										f->is_mat_known(), f->is_pat_known(),
										vector<const Progeny *>{*q});
			new_families.push_back(g);
		}
	}
	
	return new_families;
}

vector<const KnownFamily *>
				SampleManager::extract_one_known_parent_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		if((family->is_mat_known() && !is_imputed(family->get_mat()) &&
													!family->is_pat_known()) ||
				(!family->is_mat_known() && !is_imputed(family->get_pat()) &&
													family->is_pat_known())) {
			families.push_back(family);
		}
	}
	return families;
}

bool SampleManager::is_all_progenies_imputed(const KnownFamily *family) const {
	const auto&	progs = family->get_progenies();
	for(auto q = progs.begin(); q != progs.end(); ++q) {
		if(!this->is_imputed((*q)->get_name()))
			return false;
	}
	return true;
}

vector<const KnownFamily *>
			SampleManager::extract_progenies_phased_families() const {
	vector<const KnownFamily *>	families;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		if(!this->is_progeny_imputed(family))
			continue;
		
		vector<const Progeny *>	new_progenies;
		const vector<const Progeny *>&	progenies = family->get_progenies();
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			new_progenies.push_back((*q)->copy());
		}
		families.push_back(family->create(new_progenies));
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
		const KnownFamily	*family = *p;
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

vector<string> SampleManager::extract_non_imputed_samples() const {
	set<string>	samples;
	for(auto p = small_families.begin(); p != small_families.end(); ++p) {
		const KnownFamily	*family = *p;
		const auto&	progs = family->get_progenies();
		for(auto q = progs.begin(); q != progs.end(); ++q) {
			const string&	s = (*q)->get_name();
			if(s != "0" && !is_imputed(s))
				samples.insert(s);
		}
	}
	return vector<string>(samples.begin(), samples.end());
}

vector<string> SampleManager::collect_large_family_parents() const {
	set<string>	samples;
	for(auto p = large_families.begin(); p != large_families.end(); ++p) {
		const auto&	parents = (*p)->known_parents();
		for(auto q = parents.begin(); q != parents.end(); ++q)
			samples.insert(*q);
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

vector<const KnownFamily *> SampleManager::make_families(
										const PedigreeTable *ped,
										const vector<string>& samples,
										size_t lower_progs,
										const vector<size_t>& family_indices) {
	// families have been sorted
	const vector<const Family *>	families = ped->make_families(samples);
	
	vector<const KnownFamily *>	new_families;
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
		if(filtered_progs.empty())
			continue;
		
		const string&	mat = family->get_mat();
		const string&	pat = family->get_pat();
		const bool	mat_known = set_samples.find(mat) != set_samples.end();
		const bool	pat_known = set_samples.find(pat) != set_samples.end();
		if(filtered_progs.size() >= lower_progs && (mat_known || pat_known)) {
			new_families.push_back(new KnownFamily(mat, pat, mat_known,
													pat_known, filtered_progs));
		}
		else {
			// Treat parents who are not in the VCF as unknown
			const string	mat_mod = mat_known ? mat : "0";
			const string	pat_mod = pat_known ? pat : "0";
			new_families.push_back(new KnownFamily(mat_mod, pat_mod, mat_known,
													pat_known, filtered_progs));
		}
	}
	
	Common::delete_all(families);
	
	if(family_indices.empty())
		return new_families;
	
	vector<const KnownFamily *>	new_families2;
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

SampleManager *SampleManager::create(const PedigreeTable *ped,
									const vector<string>& samples,
									size_t lower_progs,
									const vector<size_t>& family_indices) {
	const auto	families = make_families(ped, samples,
											lower_progs, family_indices);
	
	vector<const KnownFamily *>	large_families;
	vector<const KnownFamily *>	small_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily *family	= *p;
		if(family->num_progenies() >= lower_progs && family->is_any_known())
			large_families.push_back(family);
		else
			small_families.push_back(family);
	}
	
	if(large_families.empty())
		throw SampleException();
	
	return new SampleManager(ped, large_families, small_families, lower_progs);
}
