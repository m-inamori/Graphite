#include <iostream>
#include <cassert>
#include "impute_VCF.h"
#include "VCFOriginal.h"
#include "VCFHeteroHomo.h"
#include "VCFCollection.h"
#include "VCFFillable.h"
#include "BiasProbability.h"
#include "option.h"
#include "common.h"

using namespace std;


//////////////////// Materials ////////////////////

Materials::~Materials() {
	delete pedigree;
	delete geno_map;
}

void Materials::select_families(set<pair<string,string>>& set_families) {
	vector<pair<string,string>>	old_families = families;
	families.clear();
	for(auto p = old_families.begin(); p != old_families.end(); ++p) {
		if(set_families.find(*p) != set_families.end())
			families.push_back(*p);
	}
}

Materials *Materials::create(const Option *option) {
	auto	*vcf = VCFOriginal::read(option->path_vcf);
	const auto	*pedigree = PedigreeTable::create(option->path_ped, vcf);
	const auto	*geno_map = Map::read(option->path_map);
	auto	families = pedigree->extract_families();
	if(option->debug)
		families.resize(3U);
	delete vcf;
	return new Materials(pedigree, geno_map, families);
}


//////////////////// process ////////////////////


BiasProbability	*bias_probability = BiasProbability::getInstance(0.01);

void select_families(Materials *materials, const HeteroParentVCFs& vcfs) {
	set<Parents>	set_families;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		set_families.insert(p->first.first);
	
	materials->select_families(set_families);
}

HeteroParentVCFs extract_VCFs(const Materials *mat, const Option *option) {
	auto	*vcf = VCFOriginal::read(option->path_vcf);
	const auto	vcfs = VCFHeteroHomo::create_vcfs(vcf, mat->get_families(),
												mat->get_ped(), mat->get_map(),
												option->debug);
	delete vcf;
	return vcfs;
}

VCFFamily *impute_each(const Parents& parents, const Map& gmap,
										const HeteroParentVCFs& vcfs, int T) {
cout << parents.first << " " << parents.second << endl;
	const int	MIN_CROSSOVER = 1;
	auto	iter_mat_vcf = vcfs.find(make_pair(parents, true));
	auto	iter_pat_vcf = vcfs.find(make_pair(parents, false));
	assert(iter_mat_vcf != vcfs.end() && iter_pat_vcf != vcfs.end());
	return VCFCollection::impute_family_vcf(iter_mat_vcf->second,
											iter_pat_vcf->second,
											gmap, MIN_CROSSOVER, T);
}

void impute_by_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	auto&	families = c->families;
	const size_t	num = families.size();
	for(size_t i = c->first; i < num; i += c->num_threads) {
		auto	*vcf = impute_each(families[i], c->gmap,
										c->vcfs, c->num_threads);
		c->imputed_vcfs[i] = vcf;
	}
}

vector<VCFFamily *> impute(const HeteroParentVCFs& vcfs,
									Materials *mat, const Option *option) {
	const auto	families = mat->get_families();
	vector<VCFFamily *>	imputed_family_vcfs(families.size());
	
	const int	T = option->num_threads;
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(families, vcfs, mat->get_map(), i, T,
														imputed_family_vcfs);
	
//#ifndef DEBUG
#if 0
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&impute_by_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_by_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	return imputed_family_vcfs;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(1);
	}
	
	Materials	*materials = Materials::create(option);
	const auto	vcfs = extract_VCFs(materials, option);
	select_families(materials, vcfs);
	
	auto	imputed_family_vcfs = impute(vcfs, materials, option);
	
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		delete p->second;
	
	const auto	*joined_vcf = VCFFillable::join_vcfs(imputed_family_vcfs,
														option->path_vcf,
														option->num_threads);
	Common::delete_all(imputed_family_vcfs);
	ofstream	ofs(option->path_out);
	joined_vcf->write(ofs);
	delete joined_vcf;
	
	delete materials;
	delete option;
}
