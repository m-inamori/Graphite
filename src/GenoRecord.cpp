#include "../include/GenoRecord.h"
#include "../include/VCF.h"
#include "../include/common.h"

using namespace std;


//////////////////// GenoRecord ////////////////////

vector<int> GenoRecord::unphased_gts() const {
	const size_t	N = num_samples();
	vector<int>	gts(N);
	for(size_t i = 0; i < N; ++i) {
		gts[i] = unphased(i);
	}
	return gts;
}

void GenoRecord::write(const VCFRecord *record, ostream& os) const {
	const auto&	orig_v = record->get_v();
	STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
	for(size_t i = 0; i < geno.size(); ++i) {
		const int	gt = geno[i];
		const string&	gt_orig = orig_v[i+9];
		const string	additional = gt_orig.substr(3, gt_orig.size() - 3);
		v.push_back(Genotype::int_to_all_gt(gt) + additional);
	}
	Common::write_tsv(v, os);
}
