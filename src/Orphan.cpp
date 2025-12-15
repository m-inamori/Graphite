#include <algorithm>
#include <cassert>
#include "../include/Orphan.h"
#include "../include/VCFOrphan.h"
#include "../include/VCFOrphanRough.h"
#include "../include/Pedigree.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// Orphan ////////////////////

// Is the computational cost sufficiently small even when using ref in HMM?
bool Orphan::is_small(const vector<vector<int>>& ref_haps,
												const OptionSmall& op) {
	const size_t	M = ref_haps[0].size();
	const size_t	NH = ref_haps.size();
	const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
	return R * M < 100000000 && R < 100000;		// 10^8 & 10^5
}

// upper NH which passes is_small
size_t Orphan::compute_upper_NH(size_t M, const OptionSmall& op) {
	for(size_t NH = 2; ; ++NH) {
		const double	R = (NH * NH * (2*NH - 1)) / op.precision_ratio;
		if(!(R * M < 100000000 && R < 100000))		// 10^8 & 10^5
			return NH - 1;
	}
	return 0;	// dummy
}

VCFGenoBase *Orphan::impute(const vector<string>& samples,
								const VCFSmall *orig_vcf,
								const vector<vector<int>>& ref_haps,
								const OptionSmall& op) {
	auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
	const size_t	N = samples.size();
	const size_t	lower_NH = 10;
	const size_t	upper_NH = 20;
	const size_t	NH = compute_upper_NH(vcf->size(), op);
	if(is_small(ref_haps, op)) {
		auto	*vcf1 = new VCFOrphan(samples, vcf->get_records(),
										ref_haps, op.map, 0.01, orig_vcf);
		vcf1->impute(op.num_threads);
		cout << N << " orphan samples have been imputed." << endl;
		vcf->clear_records();
		delete vcf;
		return vcf1;
	}
	else if(NH >= lower_NH) {
		const size_t	NH2 = min(upper_NH, NH);
		// filter the reference haplotypes
		// to those similar to the sample genotypes
		vector<vector<vector<int>>>	ref_haps_table(N);
		for(size_t i = 0; i < N; ++i) {
			const auto	gts = vcf->extract_sample_genotypes(i);
			ref_haps_table[i] =
					ReferenceHaplotype::filter_haplotypes(ref_haps, gts, NH2);
		}
		auto	*vcf2 = new VCFOrphanRough(samples, vcf->get_records(),
										ref_haps_table, op.map, 0.01, orig_vcf);
		vcf2->impute(op.num_threads);
		cout << N << " orphan samples have been imputed." << endl;
		vcf->clear_records();
		delete vcf;
		return vcf2;
	}
	else {
		delete vcf;
		return NULL;
	}
}
