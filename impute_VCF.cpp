#include <iostream>
#include <cassert>
#include "VCFOriginal.h"
#include "VCFHeteroHomo.h"
#include "VCFCollection.h"
#include "VCFFillable.h"
#include "Pedigree.h"
#include "Map.h"
#include "BiasProbability.h"
#include "option.h"
#include "common.h"

using namespace std;

BiasProbability	*bias_probability = BiasProbability::getInstance(0.01);

VCFFamily *impute_each(pair<string,string> parents, const Map& gmap,
			const map<pair<pair<string,string>,bool>,VCFHeteroHomo *>& vcfs) {
cout << parents.first << " " << parents.second << endl;
	const int	MIN_CROSSOVER = 1;
	auto	iter_mat_vcf = vcfs.find(make_pair(parents, true));
	auto	iter_pat_vcf = vcfs.find(make_pair(parents, false));
	assert(iter_mat_vcf != vcfs.end() && iter_pat_vcf != vcfs.end());
	return VCFCollection::impute_family_vcf(iter_mat_vcf->second,
											iter_pat_vcf->second,
											gmap, MIN_CROSSOVER);
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(1);
	}
	
	auto	*vcf = VCFOriginal::read(option->path_vcf);
	const auto	*pedigree = PedigreeTable::create(option->path_ped, vcf);
	auto	families = pedigree->extract_families();
	if(option->debug)
		families.resize(3U);
	const auto	*geno_map = Map::read(option->path_map);
	const auto	vcfs = VCFHeteroHomo::create_vcfs(vcf, families,
										*pedigree, *geno_map, option->debug);
	delete vcf;
	
	const auto	old_families = families;
	families.clear();
	for(auto p = old_families.begin(); p != old_families.end(); ++p) {
		if(vcfs.find(pair<pair<string,string>,bool>(*p, true)) != vcfs.end())
			families.push_back(*p);
	}
	
	// あとでマルチスレッド化したい
	vector<VCFFamily *>	imputed_family_vcfs;
	for(auto p = families.begin(); p != families.end(); ++p)
		imputed_family_vcfs.push_back(impute_each(*p, *geno_map, vcfs));
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p)
		delete p->second;
	
	const auto	joined_vcf = VCFFillable::join_vcfs(imputed_family_vcfs,
															option->path_vcf);
	Common::delete_all(imputed_family_vcfs);
	ofstream	ofs(option->path_out);
	joined_vcf->write(ofs);
	delete joined_vcf;
	
	
	delete geno_map;
	delete pedigree;
	delete option;
}
