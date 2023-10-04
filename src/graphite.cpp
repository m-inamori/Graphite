#include <iostream>
#include <cassert>
#include "../include/graphite.h"
#include "../include/SampleManager.h"
#include "../include/LargeFamily.h"
#include "../include/VCFHeteroHomoPP.h"
#include "../include/VCFOneParentPhased.h"
#include "../include/VCFProgenyPhased.h"
#include "../include/VCFIsolated.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// Materials ////////////////////

Materials::Materials(const string& path, const Map *m) : path_map(path),
										geno_map(m),
										chr_maps(Map::create_chr_maps(m)) { }

Materials::~Materials() {
	delete geno_map;
	Common::delete_all(chr_maps);
}

const Map *Materials::get_chr_map(int i) const {
	if(geno_map->is_empty())
		return chr_maps[0];
	else
		return chr_maps[i];
}

double Materials::total_cM() const {
	double	length = 0.0;
	for(auto p = chr_maps.begin(); p != chr_maps.end(); ++p)
		length += (*p)->total_cM();
	return length;
}

void Materials::display_map_info() const {
	cerr << "Genetic Map : ";
	if(geno_map->is_empty()) {
		cerr << "default map(1Mbp=1cM)." << endl;
	}
	else {
		cerr << path_map << endl;
		cerr << chr_maps.size() << " chrmosomes "
								<< total_cM() << " cM." << endl;
	}
}

Materials *Materials::create(const Option *option) {
	const auto	*geno_map = Map::read(option->path_map);
	return new Materials(option->path_map, geno_map);
}


//////////////////// process ////////////////////

VCFRecord *merge_progeny_records(vector<VCFFillable *>& vcfs,
									size_t i, const STRVEC& samples) {
	const STRVEC&	v1 = vcfs.front()->get_records()[i]->get_v();
	STRVEC	v(v1.begin(), v1.begin() + 9);
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const STRVEC&	v2 = (*p)->get_records()[i]->get_v();
		v.insert(v.end(), v2.begin() + 11, v2.end());
	}
	return new VCFRecord(v, samples);
}

VCFSmall *impute_vcf_by_parents(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				const vector<const Family *>& families,
				const Map& geno_map, int num_threads) {
	vector<VCFFillable *>	vcfs = VCFHeteroHomoPP::impute_vcfs(
														orig_vcf, merged_vcf,
														families, geno_map,
														num_threads);
	STRVEC	samples;	// 子どもだけ集める
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		VCFFillable	*vcf = *p;
		STRVEC	ss = vcf->get_samples();
		samples.insert(samples.end(), ss.begin() + 2, ss.end());
	}
	
	// samplesの所有権をnew_vcfに持たせる
	auto	new_header = orig_vcf->trim_header(samples);
	vector<VCFRecord *>	empty_records;
	auto	*new_vcf = new VCFSmall(new_header, samples, empty_records);
	const auto&	samples_ = new_vcf->get_samples();
	
	vector<VCFRecord *>	merged_records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		new_vcf->add_record(merge_progeny_records(vcfs, i, samples_));
	}
	
	Common::delete_all(vcfs);
	
	return new_vcf;
}

VCFSmall *impute_vcf_by_parent(
				const VCFSmall *orig_vcf, VCFSmall *merged_vcf,
				const vector<const Family *>& families, const Map& geno_map,
				SampleManager *sample_man, int num_threads) {
	// collect not phased parents
	STRVEC	samples;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		if(sample_man->is_imputed(family->get_mat()))
			samples.push_back(family->get_pat());
		else
			samples.push_back(family->get_mat());
	}
	
	// phase not phased parents
	VCFSmall	*parents_vcf = impute_iolated_samples(orig_vcf, merged_vcf,
														sample_man, samples,
														geno_map, num_threads);
	
	// merge vcfs
	vector<VCFSmallBase *>	vcfs{ merged_vcf, parents_vcf };
	VCFSmall	*new_merged_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
	
	// impute progenies
	VCFSmall	*new_vcf = impute_vcf_by_parents(orig_vcf, new_merged_vcf,
											families, geno_map, num_threads);
	delete new_merged_vcf;
	
	// join
	vector<VCFSmallBase *>	vcfs2{ parents_vcf, new_vcf };
	VCFSmall	*vcf = VCFSmall::join(vcfs2, orig_vcf->get_samples());
	delete parents_vcf;
	delete new_vcf;
	return vcf;
}

VCFSmall *impute_one_parent_vcf(const VCFSmall *orig_vcf,
								const VCFSmall *merged_vcf,
								const vector<const Family *>& families,
								const Map& geno_map,
								SampleManager *sample_man, int num_threads) {
	STRVEC	references = sample_man->collect_large_family_parents();
	VCFSmall	*ref_vcf = merged_vcf->extract_samples(references);
	
	vector<VCFOneParentPhased *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		const bool	is_mat_phased = sample_man->is_imputed(family->get_mat());
		auto	*vcf = VCFOneParentPhased::create(family->get_samples(),
												  is_mat_phased, merged_vcf,
												  orig_vcf, geno_map, ref_vcf);
		vcfs.push_back(vcf);
	}
	
	VCFSmall	*vcf = VCFOneParentPhased::impute_all(vcfs, num_threads);
	Common::delete_all(vcfs);
	delete ref_vcf;
	return vcf;
}

VCFSmall *impute_vcf_by_progenies(const VCFSmall *orig_vcf,
								  const VCFSmall *merged_vcf,
								  const vector<const Family *>& families,
								  const Map& geno_map,
								  SampleManager *sample_man, int num_threads) {
	STRVEC	references = sample_man->collect_large_family_parents();
	VCFSmall	*ref_vcf = merged_vcf->extract_samples(references);
	
	vector<pair<const Family *, size_t>>	progeny_imputed_families;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const Family	*family = *p;
		size_t	ppi;	// phased progeny index
		for(size_t i = 2; ; ++i) {
			if(sample_man->is_imputed(family->get_samples()[i])) {
				ppi = i;
				break;
			}
		}
		
		progeny_imputed_families.push_back(make_pair(family, ppi));
	}
	const auto	vcfs = VCFProgenyPhased::impute_all_by_progeny(
												orig_vcf, merged_vcf,
												progeny_imputed_families,
												geno_map, ref_vcf, num_threads);
	
	vector<VCFSmallBase *>	vcfs2(vcfs.begin(), vcfs.end());	// convert class
	VCFSmall	*new_vcf = VCFSmall::join(vcfs2, orig_vcf->get_samples());
	Common::delete_all(vcfs);
	delete ref_vcf;
	return new_vcf;
}

VCFSmall *impute_iolated_samples(
				const VCFSmall *orig_vcf, const VCFSmall *merged_vcf,
				SampleManager *sample_man, const STRVEC& samples,
				const Map& gmap, int num_threads) {
	const STRVEC	references = sample_man->get_large_parents();
	// あとでマルチプロセス化するためにphasingすべきsample分割する
	auto	vcfs = VCFIsolated::create(orig_vcf, merged_vcf,
										samples, references, gmap, num_threads);
	const auto	new_vcfs = VCFIsolated::impute_all(vcfs, num_threads);
	VCFSmall	*new_vcf = VCFSmall::join(new_vcfs, orig_vcf->get_samples());
	Common::delete_all(new_vcfs);
	return new_vcf;
}

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
	
	const auto&	all_samples = orig_vcf->get_samples();
	const auto	large_families = sample_man->get_large_families();
	auto	merged_vcf = LargeFamily::correct_large_family_VCFs(
									orig_vcf, large_families, geno_map, option);
	if(option->only_large_families)
		return merged_vcf;
	
	// 両親が補完されているが子どもが少ない家系を補完する
	// 補完できる家系がなくなるまで繰り返す
	sample_man->add_imputed_samples(merged_vcf->get_samples());
	while(true) {
		// 両親が補完されているが子どもが少ない家系を補完する
		auto	families1 = sample_man->extract_small_families();
		if(!families1.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_parents(orig_vcf,
												merged_vcf, families1,
												geno_map, option->num_threads);
			Common::delete_all(families1);
			vector<VCFSmallBase *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			continue;
		}
		
		// 片親が補完されている家系を補完する
		auto	families2 = sample_man->extract_single_parent_phased_families();
		if(!families2.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_parent(orig_vcf,
											merged_vcf, families2, geno_map,
											sample_man, option->num_threads);
			vector<VCFSmallBase *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			Common::delete_all(families2);
			continue;
		}
		
		auto	families3 = sample_man->extract_phased_and_unknown_parents_family();
		if(!families3.empty()) {
			auto	*new_imputed_vcf = impute_one_parent_vcf(orig_vcf,
												merged_vcf, families3, geno_map,
												sample_man, option->num_threads);
			vector<VCFSmallBase *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			Common::delete_all(families3);
			continue;
		}
		
		// Impute families whose progenies have been imputed
		auto	families4 = sample_man->extract_progenies_phased_families();
		if(!families4.empty()) {
			auto	*new_imputed_vcf = impute_vcf_by_progenies(
														orig_vcf, merged_vcf,
														families4, geno_map,
														sample_man,
														option->num_threads);
			vector<VCFSmallBase *>	vcfs{ merged_vcf, new_imputed_vcf };
			auto	*new_merged_vcf = VCFSmall::join(vcfs, all_samples);
			delete merged_vcf;
			merged_vcf = new_merged_vcf;
			sample_man->add_imputed_samples(new_imputed_vcf->get_samples());
			delete new_imputed_vcf;
			Common::delete_all(families4);
			continue;
		}
		else
			break;
	}
	
	// At last, impute isolated samples
	const STRVEC	samples = sample_man->extract_isolated_samples();
	if(!samples.empty()) {
		VCFSmall	*new_imputed_vcf = impute_iolated_samples(
												orig_vcf, merged_vcf, sample_man,
												samples,
												geno_map, option->num_threads);
		vector<VCFSmallBase *>	vcfs{ merged_vcf, new_imputed_vcf };
		VCFSmall	*vcf = merged_vcf;
		merged_vcf = VCFSmall::join(vcfs, orig_vcf->get_samples());
		delete vcf;
		delete new_imputed_vcf;
	}
	
	sample_man->clear_imputed_samples();
	return merged_vcf;
}

void print_info(const Option *option) {
	cerr << "input VCF : " << option->path_vcf << endl;
	cerr << "pedigree : " << option->path_ped << endl;
	cerr << "output VCF : " << option->path_out << endl;
}

void impute_VCF(const Option *option) {
	print_info(option);
	Materials	*materials = Materials::create(option);
	materials->display_map_info();
	
	VCFHuge	*vcf = VCFHuge::read(option->path_vcf);
	SampleManager	*sample_man = SampleManager::create(
										option->path_ped, vcf->get_samples(),
										option->lower_progs, option->families);
	sample_man->display_info();
	
	// process chromosome by chromosome
	VCFHuge::ChromDivisor	divisor(vcf);
	bool	first_chromosome = true;
	for(int	chrom_index = 0; ; ++chrom_index) {
		// chrom_index is required only for because they can skip chromosomes
		VCFSmall	*vcf_chrom = divisor.next();
		if(vcf_chrom == NULL)
			break;
		else if(!option->is_efficient_chrom(chrom_index)) {
			delete vcf_chrom;
			continue;
		}
		
		const Map	*gmap = materials->get_chr_map(chrom_index);
		const VCFSmall	*vcf_imputed = impute_vcf_chr(vcf_chrom, sample_man,
																*gmap, option);
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
	
	delete vcf;
	delete sample_man;
	delete materials;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(1);
	}
	
	impute_VCF(option);
	delete option;
}
