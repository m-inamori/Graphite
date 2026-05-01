#include <memory>

#include "../include/ImputeNoRef.h"
#include "../include/VCF.h"
#include "../include/VCFGeno.h"
#include "../include/SampleManager.h"
#include "../include/LargeFamily.h"
#include "../include/LargeSelfFamily.h"
#include "../include/materials.h"
#include "../include/SmallFamily.h"
#include "../include/option.h"

using namespace std;


//////////////////// ImputeNoRef ////////////////////

void ImputeNoRef::display_chromosome_info(const VCFSmall *orig_vcf) {
	cerr << "chr : " << orig_vcf->get_records().front()->chrom() << " ";
	if(orig_vcf->size() == 1) {
		cerr << "1 record." << endl;
	}
	else {
		cerr << orig_vcf->size() << " records." << endl;
	}
}

VCFGeno *ImputeNoRef::impute_vcf_chr(const VCFSmall *orig_vcf,
										SampleManager *sample_man,
										const Map& geno_map,
										const Option& option) {
	display_chromosome_info(orig_vcf);
	
	const auto	large_families = sample_man->get_large_families();
	auto	merged_vcf = LargeFamily::impute(orig_vcf, large_families,
															geno_map, option);
	if(merged_vcf != NULL)
		sample_man->add_imputed_samples(merged_vcf->get_samples());
	
	const auto	self_families =
						sample_man->extract_self_parent_non_imputed_families();
	auto	*imputed_vcf = LargeSelfFamily::impute(orig_vcf, merged_vcf,
											self_families, geno_map, option);
	if(imputed_vcf != NULL) {
		merged_vcf = imputed_vcf;
		sample_man->add_imputed_samples(merged_vcf->get_samples());
	}
	
	merged_vcf = SmallFamily::impute_small_family(orig_vcf, merged_vcf,
												geno_map, option, sample_man);
	
	sample_man->clear_imputed_samples();
	return merged_vcf;
}

void ImputeNoRef::impute(VCFHuge *vcf, const Materials *materials,
												const Option& option) {
	const vector<string>&	samples = vcf->get_samples();
	std::unique_ptr<const PedigreeTable> ped = nullptr;
	try {
		ped.reset(materials->get_ped()->limit_samples(samples));
	}
	catch(const ExceptionWithCode& e) {
		throw;
	}
	
	SampleManager	*sample_man = SampleManager::create(ped.get(), samples,
														vector<string>(),
														option.lower_progs,
														option.families);
	
	sample_man->display_info();
	
	// process chromosome by chromosome
	VCFHuge::ChromDivisor	divisor(vcf);
	bool	first_chromosome = true;
	for(int	chrom_index = 0; ; ++chrom_index) {
		// chrom_index is required only for because they can skip chromosomes
		const VCFSmall	*vcf_chrom = divisor.next();
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
		if(!option.is_efficient_chrom(chrom_index)) {
			delete vcf_chrom;
			continue;
		}
		
		const Map	*gmap = materials->get_chr_map(chrom_index);
		const VCFGeno	*vcf_imputed = impute_vcf_chr(vcf_chrom,
													sample_man, *gmap, option);
		if(first_chromosome) {
			ofstream	ofs(option.path_out);
			vcf_imputed->write(ofs, true);	// write header
		}
		else {
			ofstream	ofs(option.path_out, ios_base::app);
			vcf_imputed->write(ofs, false);
		}
		first_chromosome = false;
		delete vcf_imputed;
		delete vcf_chrom;
	}
	
	delete sample_man;
}
