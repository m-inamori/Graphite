#include "../include/LargeSelfFamily.h"
#include "../include/VCFSelfHetero.h"
#include "../include/VCFSelfFillable.h"
#include "../include/VCFImpSelfRecord.h"
#include "../include/VCFSelfHomoRecord.h"
#include "../include/VCFSelfJunkRecord.h"
#include "../include/SampleManager.h"
#include "../include/Pedigree.h"
#include "../include/KnownFamily.h"
#include "../include/Map.h"
#include "../include/option.h"
#include "../include/common.h"

using namespace std;


//////////////////// LargeFamily ////////////////////

// create new records and divide them into two
pair<vector<VCFSelfHeteroRecord *>, vector<VCFImpSelfRecord *>>
LargeSelfFamily::divide_records(const VCFSmall *vcf, const Option *op) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const auto	*td = CR->get_TypeDeterminer(vcf->num_samples()-1, op->ratio);
	vector<VCFSelfHeteroRecord *>	he_records;
	vector<VCFImpSelfRecord *>	other_records;
	const auto&	samples = vcf->get_samples();
	for(size_t i = 0; i < vcf->size(); ++i) {
		const VCFRecord	*record = vcf->get_record(i);
		const STRVEC&	v = record->get_v();
		const auto	pair1 = CR->classify_self_record(record, td);
		const ParentComb	pc = pair1.first;
		const WrongType	wrong_type = pair1.second;
		if(pc == ParentComb::P01x01) {
			auto	*record1 = new VCFSelfHeteroRecord(v, samples, i,
															wrong_type, pc);
			he_records.push_back(record1);
		}
		else if(TypeDeterminer::is_homohomo(pc)) {
			auto	*record2 = new VCFSelfHomoRecord(v, samples, i,
															wrong_type, pc);
			other_records.push_back(record2->impute());
			delete record2;
		}
		else {
			auto	*record3 = new VCFSelfJunkRecord(v, samples, i, wrong_type);
			other_records.push_back(record3);
		}
	}
	
	return make_pair(he_records, other_records);
}

VCFSmall *LargeSelfFamily::extract_parents(
									const vector<VCFSelfFillable *>& vcfs) {
	STRVEC	samples;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		samples.push_back((*p)->get_samples()[0]);
	}
	
	const auto	header = vcfs[0]->trim_header(samples);
	const size_t	M = vcfs[0]->size();
	vector<VCFRecord *>	records;
	for(size_t i = 0; i < M; ++i) {
		const STRVEC&	v0 = vcfs[0]->get_record(i)->get_v();
		STRVEC	v(v0.begin(), v0.begin() + 10);
		for(size_t j = 1; j < vcfs.size(); ++j) {
			v.push_back(vcfs[j]->get_record(i)->get_v()[9]);
		}
		auto	record = new VCFRecord(v, samples);
		records.push_back(record);
	}
	return new VCFSmall(header, samples, records);
}

VCFSmall *LargeSelfFamily::impute(const VCFSmall *orig_vcf,
									const VCFSmall *merged_vcf,
									const vector<const KnownFamily *>& families,
									const Map& geno_map, const Option *op) {
	if(families.empty())
		return NULL;
	
	vector<VCFSelfFillable *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		STRVEC	samples = family->get_samples();
		samples.erase(samples.begin() + 1);		// erase pat
		auto	*vcf = orig_vcf->extract_samples(samples);
		auto	pair1 = divide_records(vcf, op);
		const vector<VCFSelfHeteroRecord *>&	he_records = pair1.first;
		vector<VCFImpSelfRecord *>	other_records = pair1.second;
		auto	*vcf_hetero = new VCFSelfHetero(vcf->get_header(), samples,
														he_records, geno_map);
		
		auto	pair2 = vcf_hetero->impute(op->num_threads);
		const vector<VCFSelfHetero *>&	vcf_heteros = pair2.first;
		const vector<VCFSelfHeteroRecord *>&	unused = pair2.second;
		other_records.insert(other_records.end(), unused.begin(), unused.end());
		auto	*vcf_filled = VCFSelfFillable::fill(vcf_heteros, other_records);
		vcfs.push_back(vcf_filled);
		Common::delete_all(vcf_heteros);
		Common::delete_all(other_records);
		delete vcf_hetero;
		delete vcf;
	}
	
	const VCFSmall	*vcf_parents = extract_parents(vcfs);
	Common::delete_all(vcfs);
	auto	*new_merged_vcf = VCFSmall::join(merged_vcf, vcf_parents,
													orig_vcf->get_samples());
	delete vcf_parents;
	delete merged_vcf;
	return new_merged_vcf;
}
