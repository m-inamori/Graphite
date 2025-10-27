#include "../include/SmallFamily.h"
#include "../include/VCF.h"
#include "../include/VCFFillable.h"
#include "../include/VCFBothParentImputed.h"
#include "../include/VCFIsolated.h"
#include "../include/BothImputedFamily.h"
#include "../include/OneImputedFamily.h"
#include "../include/BothKnownFamily.h"
#include "../include/ImputedAndKnownFamily.h"
#include "../include/ProgenyImputedFamily.h"
#include "../include/OneKnownFamily.h"
#include "../include/SelfFamily.h"
#include "../include/SelfNonImputedFamily.h"
#include "../include/Orphan.h"
#include "../include/SampleManager.h"
#include "../include/KnownFamily.h"
#include "../include/Option.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;

VCFGeno *SmallFamily::impute_vcf_by_both_imputed_parents(
					const VCFSmall *orig_vcf, const VCFGenoBase *merged_vcf,
					SampleManager *sample_man, const OptionSmall& op_small) {
	auto	families = sample_man->extract_both_imputed_families();
	auto	*vcf = BothImputedFamily::impute(orig_vcf, merged_vcf,
														families, op_small);
	Common::delete_all(families);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_vcf_by_imputed_and_known_parent(
											const VCFSmall *orig_vcf,
											const VCFGeno *imputed_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	auto	families = sample_man->extract_imputed_and_known_families();
	if(families.empty())
		return NULL;
	
	// families have already been selected
	// in which one parent has been imputed and one parent has not been imputed
	// collect not phased parents
	STRVEC	samples;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*f = *p;
		if(!sample_man->is_imputed(f->get_mat()))
			samples.push_back(f->get_mat());
		if(!sample_man->is_imputed(f->get_pat()))
			samples.push_back(f->get_pat());
	}
	
	auto	*vcf = ImputedAndKnownFamily::impute_by_parent(orig_vcf,
												imputed_vcf, ref_haps,
												families, samples, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(imputed_vcf, vcf,
													orig_vcf->get_samples());
	delete imputed_vcf;
	delete vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_vcf_by_both_known_parents(
											const VCFSmall *orig_vcf,
											const VCFGeno *imputed_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	const auto	families = sample_man->extract_both_known_families();
	auto	*vcf = BothKnownFamily::impute(orig_vcf, imputed_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(imputed_vcf, vcf,
													orig_vcf->get_samples());
	delete imputed_vcf;
	delete vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_vcf_by_imputed_parent(const VCFSmall *orig_vcf,
											const VCFGeno *merged_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	const auto	families = sample_man->extract_one_imputed_families();
	auto	*vcf = OneImputedFamily::impute(orig_vcf, merged_vcf,
												ref_haps, families, op_small);
	Common::delete_all(families);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	delete vcf;
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_self_vcf(const VCFSmall *orig_vcf,
										const VCFGeno *merged_vcf,
										const vector<vector<int>>& ref_haps,
										SampleManager *sample_man,
										const OptionSmall& op_small) {
	const auto	families = sample_man->extract_self_families();
	auto	imputed_samples = sample_man->collect_imputed_samples(families);
	auto	*vcf = SelfFamily::impute(orig_vcf, merged_vcf, ref_haps,
										families, imputed_samples, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	delete vcf;
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_self_non_imputed_vcf(const VCFSmall *orig_vcf,
											const VCFGeno *merged_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	const auto	families = sample_man->extract_self_families();
	auto	*vcf = SelfNonImputedFamily::impute(orig_vcf, merged_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	delete vcf;
	return new_merged_vcf;
}

vector<vector<string>> SmallFamily::collect_imputed_progenies(
										vector<const KnownFamily *>& families,
										SampleManager *sample_man) {
	vector<vector<string>>	imputed_progenies;
	for(auto p = families.begin(); p != families.end(); ++p) {
		vector<string>	imputed_progs;
		const KnownFamily	*family = *p;
		const auto&	progs = family->get_progenies();
		for(auto q = progs.begin(); q != progs.end(); ++q) {
			const string&	prog = (*q)->get_name();
			if(sample_man->is_imputed(prog)) {
				imputed_progs.push_back(prog);
			}
		}
		if(!imputed_progs.empty())
			imputed_progenies.push_back(imputed_progs);
	}
	return imputed_progenies;
}

// Impute families whose progenies have been imputed
VCFGeno *SmallFamily::impute_vcf_by_progenies(const VCFSmall *orig_vcf,
											const VCFGeno *merged_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	auto	families = sample_man->extract_progenies_phased_families();
	if(families.empty())
		return NULL;
	
	const auto	imputed_progenies = collect_imputed_progenies(families,
																sample_man);
	auto	*vcf = ProgenyImputedFamily::impute(orig_vcf, merged_vcf,
													families, imputed_progenies,
													ref_haps, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_vcf_by_known_parent(const VCFSmall *orig_vcf,
											const VCFGeno *merged_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	auto	families = sample_man->extract_one_known_parent_families();
	auto	*vcf = OneKnownFamily::impute(orig_vcf, merged_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL)
		return NULL;
	
	auto	*new_merged_vcf = VCFGeno::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
//	Common::delete_all(families);
	return new_merged_vcf;
}

VCFGeno *SmallFamily::impute_orphan_samples(const VCFSmall *orig_vcf,
											VCFGeno *merged_vcf,
											const vector<vector<int>>& ref_haps,
											SampleManager *sample_man,
											const OptionSmall& op_small) {
	auto	samples = sample_man->extract_non_imputed_samples();
	if(samples.empty())
		return NULL;
	
	auto	*vcf = Orphan::impute(samples, orig_vcf, ref_haps, op_small);
	if(vcf == NULL)
		return NULL;
	
	// Currently, all samples may or not may be used in VCF,
	// but to be on the safe side, use samples of vcf.
	auto	*merged_vcf2 = VCFGeno::join(merged_vcf, vcf,
											orig_vcf->get_samples());
	delete merged_vcf;
	merged_vcf = merged_vcf2;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
	return merged_vcf;
}

vector<vector<int>> SmallFamily::extract_haplotypes(
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

VCFGeno *SmallFamily::impute_non_imputed_samples(const VCFSmall *orig_vcf,
													VCFGeno *merged_vcf,
													SampleManager *sample_man,
													const OptionSmall& op) {
	const STRVEC	samples = sample_man->extract_non_imputed_samples();
	if(!samples.empty() && op.imputes_isolated_samples) {
		const STRVEC	references = sample_man->collect_large_family_parents();
		// Split sample to phase for later multithreading
		auto	*vcf = VCFIsolated::create(orig_vcf, merged_vcf, samples,
																references, op);
		vcf->impute();
		VCFGeno	*new_vcf = VCFGeno::join(merged_vcf, vcf,
												orig_vcf->get_samples());
		delete merged_vcf;
		delete vcf;
		return new_vcf;
	}
	else if(!samples.empty() && op.outputs_unimputed_samples) {
		auto	*vcf_isolated = VCFGeno::extract_samples(samples, orig_vcf);
		auto	*new_vcf = VCFGeno::join(merged_vcf, vcf_isolated,
												orig_vcf->get_samples());
		delete merged_vcf;
		delete vcf_isolated;
		return new_vcf;
	}
	else {
		return merged_vcf;
	}
}

VCFGeno *SmallFamily::impute_small_family(const VCFSmall *orig_vcf,
											VCFGeno *merged_vcf,
											const Map& geno_map,
											const Option *option,
											SampleManager *sample_man) {
	OptionSmall	op_small(geno_map, option->num_threads,
										option->precision_ratio,
										option->imputes_isolated_samples,
										option->outputs_unimputed_samples);
	const auto	ref_haps = extract_haplotypes(merged_vcf, sample_man);
	// Repeat until there are no more families to impute
	while(true) {
		// Impute families in which parents are imputed but children are few
		auto	new_merged_vcf1 = impute_vcf_by_both_imputed_parents(
														orig_vcf, merged_vcf,
														sample_man, op_small);
		if(new_merged_vcf1 != NULL) {
			merged_vcf = new_merged_vcf1;
			continue;
		}
		
		// Impute families in which one parent is imputed
		auto	new_merged_vcf2 = impute_vcf_by_imputed_and_known_parent(
														orig_vcf, merged_vcf,
														ref_haps,
														sample_man, op_small);
		if(new_merged_vcf2 != NULL) {
			merged_vcf = new_merged_vcf2;
			continue;
		}
		
		auto	*new_merged_vcf3 = impute_vcf_by_both_known_parents(orig_vcf,
														merged_vcf, ref_haps,
														sample_man, op_small);
		if(new_merged_vcf3 != NULL) {
			merged_vcf = new_merged_vcf3;
			continue;
		}
		
		auto	*new_merged_vcf4 = impute_vcf_by_imputed_parent(orig_vcf,
														merged_vcf, ref_haps,
														sample_man, op_small);
		if(new_merged_vcf4 != NULL) {
			merged_vcf = new_merged_vcf4;
			continue;
		}
		
		// Impute families whose progenies have been imputed
		auto	*new_merged_vcf5 = impute_vcf_by_progenies(orig_vcf,
														merged_vcf, ref_haps,
														sample_man, op_small);
		if(new_merged_vcf5 != NULL) {
			merged_vcf = new_merged_vcf5;
			continue;
		}
		
		// Impute families in which one parent is known but not imputed
		// and the other parent is unknown
		auto	*new_merged_vcf6 = impute_vcf_by_known_parent(orig_vcf,
													merged_vcf, ref_haps,
													sample_man, op_small);
		if(new_merged_vcf6 != NULL) {
			merged_vcf = new_merged_vcf6;
			continue;
		}
		
		auto	*new_merged_vcf7 = impute_self_vcf(orig_vcf,
													merged_vcf, ref_haps,
													sample_man, op_small);
		if(new_merged_vcf7 != NULL) {
			merged_vcf = new_merged_vcf7;
			continue;
		}
		
		auto	*new_merged_vcf8 = impute_self_non_imputed_vcf(orig_vcf,
													merged_vcf, ref_haps,
													sample_man, op_small);
		if(new_merged_vcf8 != NULL) {
			merged_vcf = new_merged_vcf8;
			continue;
		}
		
		if(op_small.imputes_isolated_samples) {
			auto	*new_merged_vcf9 = impute_orphan_samples(orig_vcf,
														merged_vcf, ref_haps,
														sample_man, op_small);
			if(new_merged_vcf9 != NULL) {
				merged_vcf = new_merged_vcf9;
				continue;
			}
		}
		
		break;
	}
	
	merged_vcf = impute_non_imputed_samples(orig_vcf, merged_vcf,
													sample_man, op_small);
	return merged_vcf;
}
