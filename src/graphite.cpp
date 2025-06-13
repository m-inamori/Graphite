#include <iostream>
#include <memory>
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
	cerr << "chr : " << orig_vcf->get_records().front()->chrom() << " ";
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
	sample_man->add_imputed_samples(merged_vcf->get_samples());
	
	merged_vcf = SmallFamily::impute_small_family(orig_vcf, merged_vcf,
												geno_map, option, sample_man);
	
	sample_man->clear_imputed_samples();
	return merged_vcf;
}

void impute_all(VCFHuge *vcf, const Materials *materials,
										const Option *option) {
	const vector<string>&	samples = vcf->get_samples();
	std::unique_ptr<const PedigreeTable> ped = nullptr;
	try {
		ped.reset(materials->get_ped()->limit_samples(samples));
	}
	catch(const ExceptionWithCode& e) {
		throw;
	}
	
	SampleManager	*sample_man = SampleManager::create(ped.get(), samples,
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
