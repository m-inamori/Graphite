#include "../include/SmallFamilyRef.h"
#include "../include/BothImputedFamilyRef.h"
#include "../include/ImputedAndKnownFamilyRef.h"
#include "../include/BothKnownFamilyRef.h"
#include "../include/OneImputedFamilyRef.h"
#include "../include/OneKnownFamilyRef.h"
#include "../include/SelfFamilyRef.h"
#include "../include/SelfNonImputedFamilyRef.h"
#include "../include/ProgenyImputedFamilyRef.h"
#include "../include/OrphanRef.h"
#include "../include/VCFIsolated.h"
#include "../include/VCFGeno.h"
#include "../include/KnownFamily.h"
#include "../include/SampleManager.h"
#include "../include/OptionSmall.h"

using namespace std;


//////////////////// SmallFamilyRef ////////////////////

VCFGeno *SmallFamilyRef::merge_vcf(VCFGeno *imputed_vcf, const VCFGenoBase *vcf,
												const vector<string>& samples) {
	if(imputed_vcf == NULL) {
		return new VCFGeno(vcf->get_samples(),
							vcf->get_geno_records(), vcf->get_ref_vcf());
	}
	else {
		std::vector<const VCFGenoBase*> vcfs = { imputed_vcf, vcf };
		return VCFGeno::join(vcfs, samples);
    }
}

VCFGeno *SmallFamilyRef::impute_vcf_by_both_imputed_parents(
												const VCFSmall *orig_vcf,
												VCFGeno *phased_vcf,
												VCFGeno *imputed_vcf,
												SampleManager *sample_man,
												const OptionSmall& op_small) {
	const auto	families = sample_man->extract_both_imputed_families();
	
	const VCFGeno	*vcf = BothImputedFamilyRef::impute(orig_vcf, phased_vcf,
															families, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_vcf_by_imputed_and_known_parent(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_imputed_and_known_families();
	
	vector<string> samples;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		if(!sample_man->is_imputed(family->get_mat()))
			samples.push_back(family->get_mat());
		if(!sample_man->is_imputed(family->get_pat()))
			samples.push_back(family->get_pat());
	}
	
	VCFGeno *vcf = ImputedAndKnownFamilyRef::impute(orig_vcf, phased_vcf,
														ref_haps, families,
														samples, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_vcf_by_both_known_parents(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_both_known_families();
	VCFGeno	*vcf = BothKnownFamilyRef::impute(orig_vcf, phased_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_vcf_by_imputed_parent(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_one_imputed_families();
	
	VCFGeno *vcf = OneImputedFamilyRef::impute(orig_vcf, phased_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}


VCFGeno *SmallFamilyRef::impute_vcf_by_known_parent(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_one_known_parent_families();
	
	VCFGeno *vcf = OneKnownFamilyRef::impute(orig_vcf, phased_vcf,
												ref_haps, families, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	for(auto p = families.begin(); p != families.end(); ++p) {
		if(sample_man->is_known((*p)->get_mat()))
			sample_man->add_imputed_sample((*p)->get_mat());
		else
			sample_man->add_imputed_sample((*p)->get_pat());
	}
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_self_vcf(const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_small_self_families();
	
	vector<string> imputed_samples;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		const auto&	samples = family->get_samples();
		for(auto q = samples.begin(); q != samples.end(); ++q) {
			if(sample_man->is_imputed(*q)) {
				imputed_samples.push_back(*q);
			}
		}
	}
	
	VCFGeno *vcf = SelfFamilyRef::impute(orig_vcf, phased_vcf, ref_haps,
											families, imputed_samples, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_self_non_imputed_vcf(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_self_non_imputed_families();
	
	VCFGeno *vcf = SelfNonImputedFamilyRef::impute(orig_vcf, phased_vcf,
													ref_haps, families, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_vcf_by_progenies(
											const VCFSmall *orig_vcf,
											const VCFGeno *phased_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *imputed_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto families = sample_man->extract_progenies_phased_families();
	
	// collect imputed progenies
	vector<vector<string>> imputed_progenies;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily *family = *p;
		const auto& progenies = family->get_progenies();
		vector<string> imputed_list;
		for(auto q = progenies.begin(); q != progenies.end(); ++q) {
			if(sample_man->is_imputed((*q)->get_name())) {
				imputed_list.push_back((*q)->get_name());
			}
		}
		imputed_progenies.push_back(imputed_list);
	}
	
	VCFGeno *vcf = ProgenyImputedFamilyRef::impute(orig_vcf, phased_vcf,
													families, imputed_progenies,
													ref_haps, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(imputed_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_orphan_samples(
											const VCFSmall *orig_vcf,
											const vector<vector<int>>& ref_haps,
											VCFGeno *phased_vcf,
											const OptionSmall& op_small,
											SampleManager *sample_man) {
	const auto samples = sample_man->extract_non_imputed_samples();
	
	VCFGenoBase *vcf = OrphanRef::impute(samples, orig_vcf,
											ref_haps, phased_vcf, op_small);
	if(vcf == NULL) {
		return NULL;
	}
	
	sample_man->add_imputed_samples(vcf->get_samples());
	return merge_vcf(phased_vcf, vcf, orig_vcf->get_samples());
}

VCFGeno *SmallFamilyRef::impute_non_imputed_samples(
											const VCFSmall *orig_vcf,
											VCFGeno *merged_vcf,
											const OptionSmall& op,
											SampleManager *sample_man) {
	const auto samples = sample_man->extract_non_imputed_samples();
	
	if(!samples.empty() && !op.imputes_isolated_samples) {
		const auto	reference = sample_man->collect_reference();
		VCFIsolated *vcf = VCFIsolated::create(orig_vcf, merged_vcf,
												samples, reference, op);
		vcf->impute();
		VCFGeno *new_vcf = VCFGeno::join({merged_vcf, vcf}, orig_vcf->get_samples());
		return new_vcf;
	}
	else if(!samples.empty() && op.outputs_unimputed_samples) {
		VCFGeno *vcf_isolated = VCFGeno::extract_samples(samples, orig_vcf);
		VCFGeno *new_vcf = VCFGeno::join({merged_vcf, vcf_isolated}, orig_vcf->get_samples());
		return new_vcf;
	}
	else {
		return merged_vcf;
	}
}

VCFGeno *SmallFamilyRef::impute(const VCFSmall *orig_vcf,
								VCFGeno *merged_vcf,
								const VCFGeno *ref_vcf,
								const OptionSmall& op_small,
								SampleManager *sample_man,
								bool imputes_isolated_samples) {
	// Create reference haplotypes from the reference VCF
	const auto ref_haps = ref_vcf->create_ref_haps();
	// phased_vcf: merged_vcf combined with the reference panel (ref_vcf)
	VCFGeno	*phased_vcf = NULL;
	while(true) {
		// Merge the current merged VCF with reference VCF
		phased_vcf = VCFGeno::join(merged_vcf, ref_vcf,
													orig_vcf->get_samples());
		
		// Impute families where both parents are imputed but children are few
		if(VCFGeno *new_merged_vcf = impute_vcf_by_both_imputed_parents(
										orig_vcf, phased_vcf,
										merged_vcf, sample_man, op_small)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute families where one parent is imputed and the other is known
		if(VCFGeno *new_merged_vcf = impute_vcf_by_imputed_and_known_parent(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute families where both parents are known (not imputed)
		if(VCFGeno *new_merged_vcf = impute_vcf_by_both_known_parents(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute families where one parent is imputed
		if(VCFGeno *new_merged_vcf = impute_vcf_by_imputed_parent(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute families whose progenies have been imputed
		if(VCFGeno *new_merged_vcf = impute_vcf_by_progenies(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute families where one parent is known but not imputed and the other is unknown
		if(VCFGeno *new_merged_vcf = impute_vcf_by_known_parent(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute selfing families with at least one imputed sample
		if(VCFGeno *new_merged_vcf = impute_self_vcf(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute selfing families with no imputed samples
		if(VCFGeno *new_merged_vcf = impute_self_non_imputed_vcf(
										orig_vcf, phased_vcf,
										ref_haps, merged_vcf,
										op_small, sample_man)) {
			merged_vcf = new_merged_vcf;
			continue;
		}
		
		// Impute orphan (isolated) samples if enabled
		if(imputes_isolated_samples) {
			if(VCFGeno *new_merged_vcf = impute_orphan_samples(
											orig_vcf, ref_haps,
											phased_vcf,
											op_small, sample_man)) {
				merged_vcf = new_merged_vcf;
				continue;
			}
		}
		
		// Break the loop if no more families can be imputed
		break;
	}
	
	// Impute samples that are still non-imputed
	merged_vcf = impute_non_imputed_samples(orig_vcf, phased_vcf,
											op_small, sample_man);
	return merged_vcf;
}
