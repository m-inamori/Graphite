#include <algorithm>
#include <cassert>
#include "../include/VCFGeno.h"
#include "../include/SelfFamily.h"
#include "../include/VCFSelfParentImputed.h"
#include "../include/VCFSelfProgenyImputed.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// SelfFamily ////////////////////

VCFGeno *SelfFamily::impute(const VCFSmall *orig_vcf,
							const VCFGeno *imputed_vcf,
							const vector<vector<int>>& ref_haps,
							const vector<const KnownFamily *>& families,
							const vector<string>& imputed_samples,
							const OptionSmall& op) {
	const size_t	N = families.size();
	if(N == 0)
		return NULL;
	
	vector<const VCFGenoBase *>	vcfs;
	const set<string>	set_imputed_samples(imputed_samples.begin(),
											imputed_samples.end());
	for(size_t i = 0; i < N; ++i) {
		const KnownFamily	*family = families[i];
		// indices of imputed progenies
		vector<size_t>	prog_indices;
		const auto&	progs = family->get_progenies();
		for(size_t j = 0; j < progs.size(); ++j) {
			const string&	prog = progs[j]->get_name();
			if(set_imputed_samples.find(prog) != set_imputed_samples.end())
				prog_indices.push_back(j);
		}
		
		const string&	parent = family->get_mat();
		if(set_imputed_samples.find(parent) != set_imputed_samples.end()) {
			vector<string>	samples(1, parent);
			for(auto p = progs.begin(); p != progs.end(); ++p) {
				const string&	prog = (*p)->get_name();
				if(set_imputed_samples.find(prog) == set_imputed_samples.end())
					samples.push_back(prog);
			}
			auto	*vcf = VCFGeno::create_by_two_vcfs(imputed_vcf,
														orig_vcf, samples);
			auto	*vcf1 = new VCFSelfParentImputed(samples,
														vcf->get_records(),
														op.map, 0.01, orig_vcf);
				vcf1->impute();
			vcfs.push_back(vcf1);
			vcf->clear_records();
			delete vcf;
		}
		else if(!prog_indices.empty()) {
			vector<string>	samples(1, parent);
			samples.push_back(progs[prog_indices[0]]->get_name());
			for(auto p = progs.begin(); p != progs.end(); ++p) {
				const string&	prog = (*p)->get_name();
				if(set_imputed_samples.find(prog) == set_imputed_samples.end())
					samples.push_back(prog);
			}
			auto	*vcf = VCFGeno::create_by_two_vcfs(imputed_vcf,
														orig_vcf, samples);
			auto	*vcf1 = new VCFSelfProgenyImputed(samples, vcf->get_records(),
														ref_haps, 0, op.map, 0.01,
														orig_vcf);
			vcf1->impute();
			vcfs.push_back(vcf1);
			vcf->clear_records();
			delete vcf;
		}
	}
	
	if(vcfs.empty())
		return NULL;
	
	cout << vcfs.size() << " self families have been imputed." << endl;
	auto	*new_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	return new_vcf;
}
