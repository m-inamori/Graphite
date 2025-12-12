#include <algorithm>

#include "../include/VCFGeno.h"
#include "../include/SampleManager.h"
#include "../include/ReferenceHaplotype.h"
#include "../include/Genotype.h"

using namespace std;

vector<vector<int>> ReferenceHaplotype::extract_haplotypes(
											const VCFGeno *phased_vcf,
											const SampleManager *sample_man) {
	const STRVEC	reference = sample_man->collect_reference();
	const vector<size_t>	ref_columns = phased_vcf->extract_columns(
																reference);
	const vector<GenoRecord *>&	phased_records = phased_vcf->get_records();
	
	const size_t	NH = reference.size() * 2;
	const size_t	M = phased_records.size();
	
	// collect haplotypes
	vector<vector<int>>	gts(NH, vector<int>(M));
	for(size_t i = 0; i < M; ++i) {
		const GenoRecord	*record = phased_records[i];
		for(size_t h = 0; h < NH; ++h) {
			const size_t	c = ref_columns[h>>1];
			gts[h][i] = record->get_allele(c, h&1);
		}
	}
	
	const size_t	MIN_REF_NUM = 10;
	if(gts.size() <= MIN_REF_NUM)
		return gts;
	
	// rolling hash
	const size_t	K = std::min<size_t>(10, M);
	vector<vector<int>>	a(NH, vector<int>(M-K+1, 0));
	for(size_t h = 0; h < NH; ++h) {
		int	n = 0;
		for(size_t i = 0; i < M; ++i) {
			const int	gt = gts[h][i];
			if(i < K - 1) {
				n = (n << 1) | gt;
			}
			else {
				n = ((n << 1) & ((1 << 10) - 1)) | gt;
				a[h][i-K+1] = n;
			}
		}
	}
	
	// Check which haplotype and hash value are the same.
	vector<vector<size_t>>	b(NH, vector<size_t>(M-K+1));
	for(size_t i = 0; i < M-K+1; ++i) {
		for(size_t h = 0; h < NH; ++h) {
			b[h][i] = h;
			for(size_t l = 0; l < h; ++l) {
				if(a[l][i] == a[h][i]) {
					b[h][i] = l;
					break;
				}
			}
		}
	}
	
	// Even if the hash values are different,
	// revert to the genotype
	// and try to be as different from yourself as possible.
	for(size_t k = 1; k < NH; ++k) {
		// Group the same values
		// and convert them into data of values and ranges.
		vector<tuple<size_t, size_t, size_t>>	c;	// [(index, first, last)]
		size_t	first = 0;
		size_t	prev_h = b[k][0];
		for(size_t i = 1; i < b[k].size(); ++i) {
			const size_t	h = b[k][i];
			if(h != prev_h) {
				c.push_back(make_tuple(prev_h, first, i));
				first = i;
				prev_h = h;
			}
		}
		c.push_back(make_tuple(prev_h, first, b[k].size()));
		
		// When the haplotype itself is different before and after,
		// extend itself if possible.
		for(size_t j = 0; j < c.size(); ++j) {
			const size_t	h0     = get<0>(c[j]);
			const size_t	first0 = get<1>(c[j]);
			const size_t	last0  = get<2>(c[j]);
			if(h0 != k)
				continue;
			if(j != 0) {
				// look backward
				const size_t	h1 = get<0>(c[j-1]);
				for(size_t i = first0; i < last0; ++i) {
					if(gts[h1][i] != gts[k][i])
						break;
					b[k][i] = h1;
				}
			}
			if(j != c.size() - 1) {
				// look forward
				const size_t	h2 = get<0>(c[j+1]);
				for(size_t i = last0-1; ; --i) {
					if(gts[h2][i] != gts[k][i])
						break;
					b[k][i] = h2;
					if(i == first0)
						break;
				}
			}
		}
	}
	
	// If the difference between the two haplotypes is less than 10% in length,
	// discard one haplotype.
	// But ensure that the number does not fall below MIN_REF_NUM.
	map<size_t, size_t>	counter;
	for(size_t i = 0; i < NH; ++i) {
		for(auto p = b[i].begin(); p != b[i].end(); ++p)
			counter[*p] += 1;
	}
	
	vector<pair<size_t, size_t>>	w(NH);
	for(size_t i = 0; i < NH; ++i) {
		w[i] = make_pair(counter[i], i);
	}
	std::sort(w.begin(), w.end());
	
	vector<vector<int>>	ref_gts;
	if(w[NH-MIN_REF_NUM].first * 10 >= M) {
		for(size_t i = 0; i < NH; ++i) {
			if(w[i].first * 10 >= M)
				ref_gts.push_back(gts[w[i].second]);
		}
	}
	else {
		for(size_t i = NH-MIN_REF_NUM; i < NH; ++i) {
			ref_gts.push_back(gts[w[i].second]);
		}
	}
	return ref_gts;
}

int ReferenceHaplotype::count_same_alleles(const vector<int>& gts,
											const vector<int>& ref_hap) {
	int n = 0;
	for(size_t k = 0; k < gts.size(); ++k) {
		const int	gt = Genotype::unphased(gts[k]);
		const int	h = ref_hap[k];
		if(gt == Genotype::UN_00 && h == 0) {
			n += 1;
		}
		else if(gt == Genotype::UN_11 && h == 1) {
			n += 1;
		}
	}
	return n;
}

vector<vector<int>> ReferenceHaplotype::filter_haplotypes(
										const vector<vector<int>>& ref_haps,
										const vector<int>& gts, size_t upper) {
	const size_t	NH = ref_haps.size();
	vector<pair<int, size_t>>	ns(NH);		// [(num of same alleles, index)]
	for(size_t i = 0; i < NH; ++i) {
		ns[i] = make_pair(count_same_alleles(gts, ref_haps[i]), i);
	}
	std::sort(ns.begin(), ns.end());
	
	const size_t	first = max(NH, upper) - upper;
	vector<vector<int>>	filtered_ref_haps(NH - first);
	for(size_t i = first; i < NH; ++i) {
		filtered_ref_haps[i-first] = ref_haps[ns[i].second];
	}
	return filtered_ref_haps;
}
