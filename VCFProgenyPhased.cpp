#include <stack>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "common.h"
#include "VCFProgenyPhased.h"
#include "VCFFillable.h"
#include "Map.h"
#include "Pedigree.h"
#include "Imputer.h"

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
		auto	*record = family_records[l];
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