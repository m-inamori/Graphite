#include <memory>

#include "../include/ImputeByRef.h"
#include "../include/LargeFamilyRef.h"
#include "../include/LargeSelfFamilyRef.h"
#include "../include/SmallFamilyRef.h"
#include "../include/GenoRecord.h"
#include "../include/VCF.h"
#include "../include/VCFGeno.h"
#include "../include/SampleManager.h"
#include "../include/materials.h"
#include "../include/option.h"
#include "../include/OptionSmall.h"
#include "../include/common.h"

using namespace std;


//////////////////// ImputeByRef ////////////////////

void ImputeByRef::display_chromosome_info(const VCFSmall *orig_vcf) {
	cerr << "chr : " << orig_vcf->get_records().front()->chrom() << " ";
	if(orig_vcf->size() == 1) {
		cerr << "1 record." << endl;
	}
	else {
		cerr << orig_vcf->size() << " records." << endl;
	}
}

VCFGeno *ImputeByRef::remove_reference_samples(const VCFGeno *vcf,
											const vector<string>& ref_samples) {
	const set<string>	set_samples(ref_samples.begin(), ref_samples.end());
	
	// collect samples not included with ref_samples
	vector<string>	new_samples;
	const vector<string>&	samples = vcf->get_samples();
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(set_samples.find(*p) == set_samples.end())
			new_samples.push_back(*p);
	}
	
	// collect genotype columns
	vector<size_t>	cs;
	for(size_t c = 0; c != samples.size(); ++c) {
		if(set_samples.find(samples[c]) == set_samples.end())
			cs.push_back(c);
	}
	
	vector<GenoRecord *>	records;
	const auto&	orig_records = vcf->get_records();
	for(auto p = orig_records.begin(); p != orig_records.end(); ++p) {
		const auto	*record = *p;
		vector<int>	geno(cs.size());
		for(size_t i = 0; i < cs.size(); ++i) {
			geno[i] = record->get_geno(cs[i]);
		}
		records.push_back(new GenoRecord(record->get_pos(), geno));
	}
	
	return new VCFGeno(new_samples, records, vcf->get_ref_vcf());
}

VCFGeno *ImputeByRef::impute_vcf_chr(const VCFSmall *orig_vcf,
									 const VCFSmall *ref_vcf_,
									 SampleManager *sample_man,
									 const Map& gmap, const Option& option) {
	display_chromosome_info(orig_vcf);
	
	const VCFGeno	*ref_vcf = VCFGeno::convert(ref_vcf_);
	const auto	large_families = sample_man->get_large_families();
	VCFGeno	*merged_vcf = LargeFamilyRef::impute(large_families, orig_vcf,
														ref_vcf, gmap, option);
	if(merged_vcf != NULL)
		sample_man->add_imputed_samples(merged_vcf->get_samples());
	
	const auto	self_families =
						sample_man->extract_self_parent_non_imputed_families();
	auto	*imputed_vcf = LargeSelfFamilyRef::impute(self_families, orig_vcf,
														merged_vcf, ref_vcf,
														gmap, option);
	if(imputed_vcf != NULL) {
		merged_vcf = imputed_vcf;
		sample_man->add_imputed_samples(imputed_vcf->get_samples());
	}
	
	const OptionSmall	op_small(gmap, option.num_threads,
										option.precision_ratio,
										option.imputes_isolated_samples,
										option.outputs_unimputed_samples);
	if(merged_vcf != NULL) {
		sample_man->add_imputed_samples(merged_vcf->get_samples());
	}
	merged_vcf = SmallFamilyRef::impute(orig_vcf, merged_vcf,
											ref_vcf, op_small, sample_man,
											option.imputes_isolated_samples);
	merged_vcf->set_ref_vcf(ref_vcf_);
	sample_man->clear_imputed_samples();
	return remove_reference_samples(merged_vcf, ref_vcf->get_samples());
}

void ImputeByRef::impute(VCFHuge *vcf, const Materials *materials,
														const Option& option) {
	const auto&	samples = vcf->get_samples();
	std::unique_ptr<const PedigreeTable> ped = nullptr;
	try {
		ped.reset(materials->get_ped()->limit_samples(samples));
	}
	catch(const ExceptionWithCode& e) {
		throw;
	}
	
	VCFHuge	*ref_vcf = VCFHuge::read(option.path_ref_vcf);
	VCFHuge::ChromDivisor	divisor(vcf);
	VCFHuge::ChromDivisor	divisor_ref(ref_vcf);
	
	auto	*sample_man = SampleManager::create(ped.get(), samples,
												ref_vcf->get_samples(),
												option.lower_progs,
												option.families);
	sample_man->display_info();
	
	// process chromosome by chromosome
	bool	first_chromosome = true;
	for(int	chrom_index = 0; ; ++chrom_index) {
		// chrom_index is required only for because they can skip chromosomes
		const VCFSmall	*vcf_chrom = divisor.next();
		const VCFSmall	*ref_vcf_chrom = divisor_ref.next();
		if(vcf_chrom == NULL || ref_vcf_chrom == NULL)
			break;
		try {
			vcf_chrom->check_records();
			ref_vcf_chrom->check_records();
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
		const VCFGeno	*vcf_imputed = impute_vcf_chr(vcf_chrom, ref_vcf_chrom,
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
