#include "../include/SmallFamily.h"
#include "../include/VCF.h"
#include "../include/VCFFillable.h"
#include "../include/VCFBothParentImputed.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFIsolated.h"
#include "../include/BothImputedFamily.h"
#include "../include/OneImputedFamily.h"
#include "../include/OnePhasedFamily.h"
#include "../include/ProgenyImputedFamily.h"
#include "../include/OneKnownFamily.h"
#include "../include/Orphan.h"
#include "../include/NoPhasedFamily.h"
#include "../include/SampleManager.h"
#include "../include/KnownFamily.h"
#include "../include/common.h"

using namespace std;

VCFSmall *SmallFamily::impute_vcf_by_both_imputed_parents(
			const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
			SampleManager *sample_man,
			const Map& geno_map, const Option *option) {
	auto	families = sample_man->extract_both_imputed_families();
	if(families.empty())
		return NULL;
	
	auto	*new_imputed_vcf = BothImputedFamily::impute(orig_vcf, merged_vcf,
														families, geno_map,
														option->num_threads);
	Common::delete_all(families);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
	delete new_imputed_vcf;
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_imputed_and_known_parent(const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const vector<vector<int>>& ref_haps,
							const Map& gmap,
							SampleManager *sample_man, const Option *option) {
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
	
	auto	*vcf = OnePhasedFamily::impute_by_parent(orig_vcf, imputed_vcf,
											ref_haps, families, samples, gmap,
											option->num_threads);
	auto	*new_merged_vcf = VCFSmall::join(imputed_vcf, vcf,
													orig_vcf->get_samples());
	delete imputed_vcf;
	delete vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_both_known_parents(
							const VCFSmall *orig_vcf,
							const VCFSmall *imputed_vcf,
							const vector<vector<int>>& ref_haps,
							const Map& gmap,
							SampleManager *sample_man, int num_threads) {
	const auto	families = sample_man->extract_both_known_families();
	if(families.empty())
		return NULL;
	
	auto	*vcf = NoPhasedFamily::impute(orig_vcf, imputed_vcf, ref_haps,
												families, gmap, num_threads);
	auto	*new_merged_vcf = VCFSmall::join(imputed_vcf, vcf,
													orig_vcf->get_samples());
	delete imputed_vcf;
	delete vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_imputed_parent(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const vector<vector<int>>& ref_haps,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads) {
	const auto	families = sample_man->extract_one_imputed_families();
	if(families.empty())
		return NULL;
	
	auto	*vcf = OneImputedFamily::impute(orig_vcf, merged_vcf,
												ref_haps, families,
												geno_map, num_threads);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(new_merged_vcf->get_samples());
	delete vcf;
	Common::delete_all(families);
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
VCFSmall *SmallFamily::impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const vector<vector<int>>& ref_haps,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads) {
	auto	families = sample_man->extract_progenies_phased_families();
	if(families.empty())
		return NULL;
	
	const auto	imputed_progenies = collect_imputed_progenies(families,
																sample_man);
	auto	*vcf = ProgenyImputedFamily::impute(orig_vcf, merged_vcf,
													families, imputed_progenies,
													ref_haps, geno_map,
													num_threads);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_vcf_by_known_parent(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const vector<vector<int>>& ref_haps,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads) {
	auto	families = sample_man->extract_one_known_parent_families();
	if(families.empty())
		return NULL;
	
	auto	*vcf = OneKnownFamily::impute(orig_vcf, merged_vcf, ref_haps,
											families, geno_map, num_threads);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, vcf,
													orig_vcf->get_samples());
	delete merged_vcf;
	sample_man->add_imputed_samples(vcf->get_samples());
	delete vcf;
//	Common::delete_all(families);
	return new_merged_vcf;
}

VCFSmall *SmallFamily::impute_orphan_samples(const VCFSmall *orig_vcf,
								  VCFSmall *merged_vcf,
								  const vector<vector<int>>& ref_haps,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads) {
	auto	samples = sample_man->extract_non_imputed_samples();
	if(samples.empty())
		return NULL;
	
	auto	*vcf = Orphan::impute(samples, orig_vcf, ref_haps,
												geno_map, num_threads);
	if(vcf != NULL) {
		// Currently, all samples may or not may be used in VCF,
		// but to be on the safe side, use samples of vcf.
		auto	*merged_vcf2 = VCFSmall::join(merged_vcf, vcf,
												orig_vcf->get_samples());
		delete merged_vcf;
		merged_vcf = merged_vcf2;
		sample_man->add_imputed_samples(vcf->get_samples());
		delete vcf;
	}
	return merged_vcf;
}

vector<vector<int>> SmallFamily::extract_haplotypes(
										const VCFSmall *phased_vcf,
										const SampleManager *sample_man) {
	const STRVEC	reference = sample_man->collect_reference();
	const vector<size_t>	ref_columns = phased_vcf->extract_columns(
																reference);
	const vector<VCFRecord *>&	phased_records = phased_vcf->get_records();
	
	const size_t	NH = reference.size() * 2;
	const size_t	M = phased_records.size();
	
	// collect haplotypes
	vector<vector<int>>	gts(NH, vector<int>(M));
	for(size_t i = 0; i < M; ++i) {
		const STRVEC&	v = phased_records[i]->get_v();
		for(size_t k = 0; k < NH; ++k) {
			const size_t	c = ref_columns[k>>1];
			gts[k][i] = (int)(v[c].c_str()[(k&1)<<1] - '0');
		}
	}
	
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
	map<size_t, size_t>	counter;
	for(size_t i = 0; i < NH; ++i) {
		for(auto p = b[i].begin(); p != b[i].end(); ++p)
			counter[*p] += 1;
	}
	
	vector<vector<int>>	ref_gts;
	for(auto p = counter.begin(); p != counter.end(); ++p) {
		if(p->second * 10 >= M)
			ref_gts.push_back(gts[p->first]);
	}
	return ref_gts;
}

VCFSmall *SmallFamily::impute_small_family_VCFs(const VCFSmall *orig_vcf,
												VCFSmall *merged_vcf,
												const Map& geno_map,
												SampleManager *sample_man,
												const Option *option) {
	const auto	ref_haps = extract_haplotypes(merged_vcf, sample_man);
	// Repeat until there are no more families to impute
	while(true) {
		// Impute families in which parents are imputed but children are few
		auto	new_merged_vcf1 = impute_vcf_by_both_imputed_parents(orig_vcf,
														merged_vcf, sample_man,
														geno_map, option);
		if(new_merged_vcf1 != NULL) {
			merged_vcf = new_merged_vcf1;
			continue;
		}
		
		// Impute families in which one parent is imputed
		auto	new_merged_vcf2 = impute_vcf_by_imputed_and_known_parent(
														orig_vcf, merged_vcf,
														ref_haps, geno_map,
														sample_man, option);
		if(new_merged_vcf2 != NULL) {
			merged_vcf = new_merged_vcf2;
			continue;
		}
		
		auto	*new_merged_vcf3 = impute_vcf_by_both_known_parents(
														orig_vcf, merged_vcf,
														ref_haps,
														geno_map, sample_man,
														option->num_threads);
		if(new_merged_vcf3 != NULL) {
			merged_vcf = new_merged_vcf3;
			continue;
		}
		
		auto	*new_merged_vcf4 = impute_vcf_by_imputed_parent(
														orig_vcf, merged_vcf,
														ref_haps,
														geno_map, sample_man,
														option->num_threads);
		if(new_merged_vcf4 != NULL) {
			merged_vcf = new_merged_vcf4;
			continue;
		}
		
		// Impute families whose progenies have been imputed
		auto	*new_merged_vcf5 = impute_vcf_by_progenies(orig_vcf,
													merged_vcf, ref_haps,
													geno_map, sample_man,
													option->num_threads);
		if(new_merged_vcf5 != NULL) {
			merged_vcf = new_merged_vcf5;
			continue;
		}
		
		// Impute families in which one parent is known but not imputed
		// and the other parent is unknown
		auto	*new_merged_vcf6 = impute_vcf_by_known_parent(orig_vcf,
													merged_vcf, ref_haps,
													geno_map, sample_man,
													option->num_threads);
		if(new_merged_vcf6 != NULL) {
			merged_vcf = new_merged_vcf6;
			continue;
		}
		
		auto	*new_merged_vcf7 = impute_orphan_samples(orig_vcf,
													merged_vcf, ref_haps,
													geno_map, sample_man,
													option->num_threads);
		if(new_merged_vcf7 != NULL) {
			merged_vcf = new_merged_vcf7;
			continue;
		}
		
		break;
	}
	return merged_vcf;
}

VCFSmall *SmallFamily::impute_iolated_samples(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				SampleManager *sample_man, const STRVEC& samples,
				const Map& gmap, bool modify_genotypes, int num_threads) {
	const STRVEC	references = sample_man->collect_large_family_parents();
	// Split sample to phase for later multithreading
	auto	vcfs_ = VCFIsolated::create(orig_vcf, merged_vcf, samples,
										references, gmap,
										modify_genotypes, num_threads);
	const auto	new_vcfs = VCFIsolated::impute_all(vcfs_, num_threads);
	vector<const VCFSmallBase *>	vcfs(new_vcfs.begin(), new_vcfs.end());
	VCFSmall	*new_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	Common::delete_all(new_vcfs);
	return new_vcf;
}
