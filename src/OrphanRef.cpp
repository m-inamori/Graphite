#include "../include/OrphanRef.h"
#include "../include/Orphan.h"
#include "../include/VCFGeno.h"
#include "../include/RefCommon.h"

using namespace std;


//////////////////// OrphanRef ////////////////////

VCFGenoBase *OrphanRef::impute(const std::vector<std::string>& samples,
								const VCFSmall *orig_vcf,
								const std::vector<std::vector<int>>& ref_haps,
								const VCFGeno *phased_vcf,
								const OptionSmall& op) {
	if(samples.empty())
		return NULL;
	
	auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
	const auto	records = RefCommon::expand_records(vcf, phased_vcf);
	auto	*vcf1 = Orphan::impute_samples(samples, records,
											ref_haps, orig_vcf, op);
	if(vcf1 != NULL)
		cerr << samples.size() << " orphan samples have been imputed." << endl;
	delete vcf;
	return vcf1;
}
