#include <algorithm>
#include <cassert>
#include "../include/ProgenyImputedFamilyRef.h"
#include "../include/ProgenyImputedFamily.h"
#include "../include/VCFFamily.h"
#include "../include/VCFProgenyImputed.h"
#include "../include/KnownFamily.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// ProgenyImputedFamilyRef ////////////////////

VCFGeno *ProgenyImputedFamilyRef::impute(const VCFSmall *orig_vcf,
								const VCFGenoBase *imputed_vcf,
								const vector<const KnownFamily *>& families,
								const vector<vector<string>>& imputed_progenies,
								const vector<vector<int>>& ref_haps,
								const OptionSmall& op) {
	const auto	sample_table = ProgenyImputedFamily::collect_samples(families,
															imputed_progenies);
	
	vector<VCFImputable *>	vcfs;
	vector<VCFFamily *>	vcf_garbage;
	const size_t	L = families.size();
	for(size_t i = 0; i < L; ++i) {
		const KnownFamily	*family = families[i];
		auto	*vcf = VCFGeno::extract_samples(sample_table[i], orig_vcf);
		const auto	records = RefCommon::merge_family_records(imputed_vcf, vcf,
															sample_table[i]);
		auto	*vcf1 = new VCFProgenyImputed(sample_table[i], records,
											  ref_haps, family->is_mat_known(),
											  op.map, 0.01, orig_vcf);
		vcfs.push_back(vcf1);
		vcf->clear_records();
		delete vcf;
	}
	
	if(vcfs.empty()) {
		Common::delete_all(vcf_garbage);
		return NULL;
	}
	
	// Small VCFs are heavy to process, so it will be parallelized.
	VCFImputable::impute_VCFs(vcfs, op.num_threads);
	cout << vcfs.size()
			<< " families whose progeny is imputed have been imputed." << endl;
	
	vector<const VCFGenoBase *>	vcfs1(vcfs.begin(), vcfs.end());
	auto	*new_vcf = VCFGeno::join(vcfs1, orig_vcf->get_samples());
	Common::delete_all(vcf_garbage);
	Common::delete_all(vcfs);
	return new_vcf;
}
