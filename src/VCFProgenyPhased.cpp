#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFFillable.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"

using namespace std;


//////////////////// VCFProgenyPhased ////////////////////

VCFProgenyPhased::VCFProgenyPhased(const vector<STRVEC>& h, const STRVEC& s,
					vector<VCFFamilyRecord *> rs, const vector<size_t>& ppi) :
			VCFFamily(h, s, vector<VCFFamilyRecord *>(rs.begin(), rs.end())),
			phased_progeny_indices(ppi) { }

void VCFProgenyPhased::determine_parent(bool is_mat) {
	// First sample is used for now
	// I'd like to use all the samples if possible
	const size_t	i = phased_progeny_indices[0];
	const size_t	j = is_mat ? 0 : 2;
	const size_t	k = is_mat ? 0 : 1;
	for(size_t l = 0; l < this->size(); ++l) {
		auto	*record = records[l];
		const string	prog_GT = record->get_GT(i);
		const int	int_gt = record->get_int_gt(k);
		const int	remain = int_gt - (int)(prog_GT.c_str()[j] - '0');
		// Let the left transfer to the progeny
		if(remain <= 0)
			record->set_GT(k, prog_GT.substr(j, 1) + "|0");
		else
			record->set_GT(k, prog_GT.substr(j, 1) + "|1");
	}
}

void VCFProgenyPhased::impute() {
	if(this->size() == 0)
		return;
	
	// Impute only parents
	// Remaining progenies are phased in impute_vcf_by_parents
	this->determine_parent(true);
	this->determine_parent(false);
}

VCFProgenyPhased *VCFProgenyPhased::impute_by_progeny(const VCFSmall *orig_vcf,
													const VCFSmall *imputed_vcf,
													const STRVEC& samples,
													const vector<size_t>& ppi) {
	VCFFamily	*vcf = VCFFamily::create_by_two_vcfs(imputed_vcf,
														orig_vcf, samples);
	auto	new_vcf = new VCFProgenyPhased(vcf->get_header(),
							vcf->get_samples(), vcf->get_family_records(), ppi);
	new_vcf->impute();
	vcf->clear_records();
	delete vcf;
	return new_vcf;
}

void VCFProgenyPhased::impute_in_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	for(size_t i = c->first; i < c->size(); i += c->num_threads) {
		const Family	*family = c->families[i].first;
		const PPI		ppi = c->families[i].second;
		c->results[i] = impute_by_progeny(c->orig_vcf, c->merged_vcf,
												family->get_samples(), ppi);
	}
}

vector<VCFProgenyPhased *> VCFProgenyPhased::impute_all_by_progeny(
							const VCFSmall *orig_vcf,
							const VCFSmall *merged_vcf,
							const vector<pair<const Family *, PPI>>& families,
							int num_threads) {
	vector<VCFProgenyPhased *>	results(families.size());
	
	const int	T = min((int)families.size(), num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(orig_vcf, merged_vcf, families,
													(size_t)i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
						(void *(*)(void *))&impute_in_thread,
						(void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	return results;
}
