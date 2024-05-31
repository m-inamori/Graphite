#include <iostream>
#include <cassert>
#include "../include/graphite.h"
#include "../include/SampleManager.h"
#include "../include/LargeFamily.h"
#include "../include/SmallFamily.h"
#include "../include/impute_prog_only.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFIsolated.h"
#include "../include/materials.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// process ////////////////////

void display_chromosome_info(const VCFSmall *orig_vcf) {
	cerr << "chr : " << orig_vcf->get_records().front()->chrom() << endl;
	if(orig_vcf->size() == 1) {
		cerr << "1 record." << endl;
	}
	else {
		cerr << orig_vcf->size() << " records." << endl;
	}
}

VCFSmall *impute_vcf_chr(const VCFSmall *orig_vcf, SampleManager *sample_man,
									const Map& geno_map, const Option *option) {
	display_chromosome_info(orig_vcf);
	
	const auto	large_families = sample_man->get_large_families();
	auto	merged_vcf = LargeFamily::correct_large_family_VCFs(
									orig_vcf, large_families, geno_map, option);
	if(option->only_large_families)
		return merged_vcf;
	sample_man->add_imputed_samples(merged_vcf->get_samples());
	
	merged_vcf = SmallFamily::impute_small_family_VCFs(orig_vcf, merged_vcf,
														geno_map, sample_man,
														option);
	
	// At last, impute isolated samples
	const STRVEC	samples = sample_man->extract_isolated_samples();
	if(!samples.empty() && option->imputes_isolated_samples) {
		VCFSmall	*new_imputed_vcf = SmallFamily::impute_iolated_samples(
											orig_vcf, merged_vcf,
											sample_man, samples, geno_map,
											option->corrects_isolated_samples,
											option->num_threads);
		VCFSmall	*vcf = merged_vcf;
		merged_vcf = VCFSmall::join(merged_vcf, new_imputed_vcf,
												orig_vcf->get_samples());
		delete vcf;
		delete new_imputed_vcf;
	}
	else if(!samples.empty() && option->outputs_unimputed_samples) {
		auto	*vcf_isolated = orig_vcf->extract_samples(samples);
		VCFSmall	*vcf = merged_vcf;
		merged_vcf = VCFSmall::join(merged_vcf, vcf_isolated,
												orig_vcf->get_samples());
		delete vcf;
	}
	
	sample_man->clear_imputed_samples();
	return merged_vcf;
}

void impute_all(VCFHuge *vcf, const Materials *materials,
										const Option *option) {
	SampleManager	*sample_man = SampleManager::create(
									materials->get_ped(), vcf->get_samples(),
									option->lower_progs, option->families);
	
	sample_man->display_info();
	
	// process chromosome by chromosome
	VCFHuge::ChromDivisor	divisor(vcf);
	bool	first_chromosome = true;
	for(int	chrom_index = 0; ; ++chrom_index) {
		// chrom_index is required only for because they can skip chromosomes
		VCFSmall	*vcf_chrom;
		vcf_chrom = divisor.next();
		if(vcf_chrom == NULL)
			break;
		try {
			vcf_chrom->check_records();
		}
		catch(const std::exception& e) {
			delete vcf_chrom;
			delete sample_man;
			throw;
		}
		if(!option->is_efficient_chrom(chrom_index)) {
			delete vcf_chrom;
			continue;
		}
		
		const Map	*gmap = materials->get_chr_map(chrom_index);
		const VCFSmall	*vcf_imputed = impute_vcf_chr(vcf_chrom,
													sample_man, *gmap, option);
		delete vcf_chrom;
		if(first_chromosome) {
			ofstream	ofs(option->path_out);
			vcf_imputed->write(ofs, true);	// write header
		}
		else {
			ofstream	ofs(option->path_out, ios_base::app);
			vcf_imputed->write(ofs, false);
		}
		first_chromosome = false;
		delete vcf_imputed;
	}
	
	delete sample_man;
}

void impute_progenies(VCFHuge *vcf, const Materials *materials,
												const Option *option) {
	auto	*vcf_ref = VCFHuge::read(option->path_ref_vcf);
	VCFHuge::ChromDivisor	divisor(vcf);
	VCFHuge::ChromDivisor	divisor_ref(vcf_ref);
	bool	first_chromosome = true;
	// とりあえず、後代のVCFも同じ染色体があるとする
	for(int	chrom_index = 0; ; ++chrom_index) {
		// chrom_index is required only for because they can skip chromosomes
		VCFSmall	*vcf_chrom = divisor.next();
		VCFSmall	*vcf_ref_chrom = divisor_ref.next();
		if(vcf_chrom == NULL)
			break;
		else if(!option->is_efficient_chrom(chrom_index)) {
			delete vcf_chrom;
			continue;
		}
		
		const Map	*gmap = materials->get_chr_map(chrom_index);
		const VCFSmallBase	*vcf_imputed = ImputeProgOnly::impute_prog_vcf_chr(
													vcf_ref_chrom,
													vcf_chrom, *gmap, option);
		delete vcf_chrom;
		delete vcf_ref_chrom;
		if(first_chromosome) {
			ofstream	ofs(option->path_out);
			vcf_imputed->write(ofs, true);	// write header
		}
		else {
			ofstream	ofs(option->path_out, ios_base::app);
			vcf_imputed->write(ofs, false);
		}
		first_chromosome = false;
		delete vcf_imputed;
	}
}

void impute_VCF(const Option *option) {
	option->print_info();
	Materials	*materials = Materials::create(option->path_map,
												option->path_ped);
	materials->display_map_info();
	
	VCFHuge	*vcf = NULL;
	try {
		vcf = VCFHuge::read(option->path_vcf);
		if(option->exists_ref())
			impute_progenies(vcf, materials, option);
		else
			impute_all(vcf, materials, option);
		delete vcf;
	}
	catch(const ExceptionWithCode& e) {
		if(vcf != NULL)
			delete vcf;
		delete materials;
		throw;
	}
	delete materials;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(ErrorCode::WRONG_ARGUMENT);
	}
	
	try {
		impute_VCF(option);
	}
	catch(const ExceptionWithCode& e) {
		delete option;
		cerr << e.what() << endl;
		exit(e.get_error_code());
	}
	
	delete option;
	return 0;
}
