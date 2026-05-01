#include "../include/VCF.h"
#include "../include/GenoRecord.h"
#include "../include/VCFGeno.h"
#include "../include/common.h"

using namespace std;


//////////////////// GenoRecord ////////////////////

GenoRecord *GenoRecord::copy() const {
	vector<int>	new_geno(geno.size());
	std::copy(geno.begin(), geno.end(), new_geno.begin());
	return new GenoRecord(pos, new_geno);
}

vector<int> GenoRecord::unphased_gts() const {
	const size_t	N = num_samples();
	vector<int>	gts(N);
	for(size_t i = 0; i < N; ++i) {
		gts[i] = unphased(i);
	}
	return gts;
}

void GenoRecord::copy_genotypes_from(const GenoRecord *other) {
	std::copy(other->geno.begin(), other->geno.end(), geno.begin());
}

void GenoRecord::write(const VCFRecord *record,
						const vector<size_t>& columns, ostream& os) const {
	const auto&	orig_v = record->get_v();
	STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
	const STRVEC	def_info = record->default_info();
	for(size_t i = 0; i < geno.size(); ++i) {
		const int	gt = geno[i];
		// Try to retrieve the original extra information (non-GT fields) for this sample.
		// If not available, fall back to the default information prepared above.
		const STRVEC	ext_info = record->extra_info(columns[i], def_info);
		vector<string>	w(1, Genotype::int_to_all_gt(gt));
		w.insert(w.end(), ext_info.begin(), ext_info.end());
		v.push_back(Common::join(w, ':'));
	}
	Common::write_tsv(v, os);
}

vector<int> GenoRecord::extract_sample_genotypes(size_t c,
										const vector<GenoRecord *>& records) {
	vector<int>	gts(records.size());
	for(size_t i = 0; i < records.size(); ++i) {
		gts[i] = records[i]->geno[c];
	}
	return gts;
}
