#ifndef __IMPUTEVCF
#define __IMPUTEVCF

class Option;


//////////////////// process ////////////////////

// Perform genotype imputation for the input VCF.
// The imputation method is selected depending on whether
// a reference panel is available.
void graphite(const Option& option);
#endif
