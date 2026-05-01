#include "../include/RefCommon.h"
#include "../include/GenoRecord.h"
#include "../include/VCFGeno.h"
#include "../include/VCFFamily.h"

using namespace std;


//////////////////// RefCommon ////////////////////

vector<pair<size_t, vector<int>>>
RefCommon::merge_records_core(const VCFGenoBase *phased_vcf,
								const VCFGenoBase *non_phased_vcf,
								const STRVEC& samples) {
	// fill with N/A
	const size_t	N = samples.size();
	const size_t	M = phased_vcf->size();
	vector<pair<size_t, vector<int>>>	genos(M);
	for(size_t i = 0; i < M; ++i) {
		const GenoRecord	*record = phased_vcf->get_record(i);
		const ll	pos = record->get_pos();
		const vector<int>	geno(N, Genotype::NA);
		genos[i] = make_pair(pos, geno);
	}
	
	// Assign genotypes from the phased VCF
	const auto	cs1 = phased_vcf->extract_columns(samples);
	const auto	cs2 = non_phased_vcf->extract_columns(samples);
	for(size_t j = 0; j < M; ++j) {
		const GenoRecord	*record1 = phased_vcf->get_record(j);
		const auto&	geno = record1->get_genos();
		for(size_t i = 0; i < cs1.size(); ++i) {
			const size_t	c = cs1[i];
			if(c != string::npos)
				genos[j].second[i] = geno[c];
		}
	}
	
	size_t	k = 0;
	size_t	l = 0;
	const size_t	L1 = phased_vcf->size();
	const size_t	L2 = non_phased_vcf->size();
	while(k < L1 && l < L2) {
		const GenoRecord	*record1 = phased_vcf->get_record(k);
		const GenoRecord	*record2 = non_phased_vcf->get_record(l);
		if(record1->get_pos() == record2->get_pos()) {
			const auto&	geno = record2->get_genos();
			for(size_t i = 0; i < cs2.size(); ++i) {
				const size_t	c = cs2[i];
				if(c != string::npos && cs1[i] == string::npos)
					genos[k].second[i] = geno[c];
			}
		}
		if(record1->get_pos() <= record2->get_pos()) {
			k += 1;
		}
		if(record1->get_pos() >= record2->get_pos()) {
			l += 1;
		}
	}
	
	return genos;
}

vector<GenoRecord *> RefCommon::merge_records(const VCFGenoBase *phased_vcf,
											const VCFGenoBase *non_phased_vcf,
											const STRVEC& samples) {
	const auto	genos = merge_records_core(phased_vcf, non_phased_vcf, samples);
	const size_t	M = genos.size();
	vector<GenoRecord *>	records(M);
	for(size_t i = 0; i < M; ++i) {
		records[i] = new GenoRecord(genos[i].first, genos[i].second);
	}
	return records;
}

vector<VCFFamilyRecord *> RefCommon::merge_family_records(
											const VCFGenoBase *phased_vcf,
											const VCFGenoBase *non_phased_vcf,
											const STRVEC& samples) {
	const auto	genos = merge_records_core(phased_vcf, non_phased_vcf, samples);
	const size_t	M = genos.size();
	vector<VCFFamilyRecord *>	records(M);
	for(size_t i = 0; i < M; ++i) {
		records[i] = new VCFFamilyRecord(genos[i].first, genos[i].second);
	}
	return records;
}

vector<GenoRecord *> RefCommon::expand_records(const VCFGeno *vcf,
												const VCFGeno *phased_vcf) {
	if(vcf->size() == 0)
		return vector<GenoRecord *>();
	
	// fill with N/A
	const size_t	N = vcf->get_samples().size();
	const size_t	M = phased_vcf->size();
	vector<GenoRecord *>	records(M);
	for(size_t i = 0; i < M; ++i) {
		const GenoRecord	*record = phased_vcf->get_record(i);
		const ll	pos = record->get_pos();
		const vector<int>	geno(N, Genotype::NA);
		records[i] = new GenoRecord(pos, geno);
	}
	
	size_t	k = 0;
	size_t	l = 0;
	const size_t	L1 = vcf->size();
	const size_t	L2 = phased_vcf->size();
	while(k < L1 && l < L2) {
		const GenoRecord	*record1 = vcf->get_record(k);
		GenoRecord	*record2 = records[l];
		if(record1->get_pos() == record2->get_pos()) {
			record2->copy_genotypes_from(record1);
		}
		if(record1->get_pos() <= record2->get_pos()) {
			k += 1;
		}
		if(record1->get_pos() >= record2->get_pos()) {
			l += 1;
		}
	}
	
	return records;
}
