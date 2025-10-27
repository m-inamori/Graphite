#ifndef __IMPUTEVCF
#define __IMPUTEVCF

#include <map>
#include <set>
#include "VCFHeteroHomo.h"
#include "VCFFillable.h"
#include "Pedigree.h"
#include "Map.h"
#include "error_codes.h"

class Option;
class SampleManager;
class VCFFamily;
class VCFHeteroHomo;
class Materials;


//////////////////// process ////////////////////

void display_chromosome_info(const VCFSmall *orig_vcf);
VCFGeno*impute_vcf_chr(const VCFSmall *vcf, SampleManager *sample_man,
								const Map& geno_map, const Option *option);
void print_info(const Option *option);
void impute_all(VCFHuge *vcf, const Materials *materials, const Option *option);
void impute_progenies(VCFHuge *vcf, const Materials *materials,
												const Option *option);
void impute_VCF(const Option *option);
#endif
