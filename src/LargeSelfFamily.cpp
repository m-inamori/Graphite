#include "../include/VCFGeno.h"
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


//////////////////// LargeSelfFamily ////////////////////

// create new records and divide them into two
pair<vector<VCFSelfHeteroRecord *>, vector<VCFImpSelfRecord *>>
LargeSelfFamily::divide_records(const VCFGeno *vcf, const Option *op) {
	ClassifyRecord	*CR = ClassifyRecord::get_instance();
	const auto	*td = CR->get_TypeDeterminer(vcf->num_samples()-1, op->ratio);
	vector<VCFSelfHeteroRecord *>	he_records;
	vector<VCFImpSelfRecord *>	other_records;
	for(size_t i = 0; i < vcf->size(); ++i) {
		const GenoRecord	*record = vcf->get_record(i);
		const ll	pos = record->get_pos();
		const auto	geno = record->get_geno();
		const auto	pair1 = CR->classify_self_record(record, td);
		const ParentComb	pc = pair1.first;
		const WrongType	wrong_type = pair1.second;
		if(pc == ParentComb::P01x01) {
			auto	*record1 = new VCFSelfHeteroRecord(pos, geno, i,
															wrong_type, pc);
			he_records.push_back(record1);
		}
		else if(TypeDeterminer::is_homohomo(pc)) {
			auto	*record2 = new VCFSelfHomoRecord(pos, geno, i,
															wrong_type, pc);
			other_records.push_back(record2->impute());
			delete record2;
		}
		else {
			auto	*record3 = new VCFSelfJunkRecord(pos, geno, i, wrong_type);
			other_records.push_back(record3);
		}
	}
	
	return make_pair(he_records, other_records);
}

VCFGeno *LargeSelfFamily::extract_parents(
									const vector<VCFSelfFillable *>& vcfs) {
	STRVEC	samples;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		samples.push_back((*p)->get_samples()[0]);
	}
	
	const size_t	M = vcfs[0]->size();
	vector<GenoRecord *>	records;
	for(size_t i = 0; i < M; ++i) {
		const ll	pos = vcfs[0]->get_record(i)->get_pos();
		vector<int>	geno;
		for(size_t j = 0; j < vcfs.size(); ++j) {
			geno.push_back(vcfs[j]->get_record(i)->get_geno()[0]);
		}
		auto	record = new GenoRecord(pos, geno);
		records.push_back(record);
	}
	return new VCFGeno(samples, records, vcfs[0]->get_ref_vcf());
}

VCFGeno *LargeSelfFamily::impute(const VCFSmall *orig_vcf,
									VCFGeno *merged_vcf,
									const vector<const KnownFamily *>& families,
									const Map& geno_map, const Option *op) {
	if(families.empty())
		return NULL;
	
	vector<VCFSelfFillable *>	vcfs;
	for(auto p = families.begin(); p != families.end(); ++p) {
		const KnownFamily	*family = *p;
		STRVEC	samples = family->get_samples();
		samples.erase(samples.begin() + 1);		// erase pat
		auto	*vcf = VCFGeno::extract_samples(samples, orig_vcf);
		auto	pair1 = divide_records(vcf, op);
		const vector<VCFSelfHeteroRecord *>&	he_records = pair1.first;
		vector<VCFImpSelfRecord *>	other_records = pair1.second;
		auto	*vcf_hetero = new VCFSelfHetero(samples, he_records,
													geno_map, 0.01, orig_vcf);
		
		auto	pair2 = vcf_hetero->impute(op->num_threads);
		const vector<VCFSelfHetero *>&	vcf_heteros = pair2.first;
		const vector<VCFSelfHeteroRecord *>&	unused = pair2.second;
		other_records.insert(other_records.end(), unused.begin(), unused.end());
		auto	*vcf_filled = VCFSelfFillable::fill(vcf_heteros, other_records);
		vcfs.push_back(vcf_filled);
		Common::delete_all(vcf_heteros);
		Common::delete_all(other_records);
		vcf_hetero->clear_records();
		delete vcf_hetero;
		delete vcf;
	}
	
	const VCFGeno	*vcf_parents = extract_parents(vcfs);
	Common::delete_all(vcfs);
	VCFGeno	*new_merged_vcf;
	if(merged_vcf != NULL) {
		new_merged_vcf = VCFGeno::join(merged_vcf, vcf_parents,
													orig_vcf->get_samples());
	}
	else {
		vector<const VCFGenoBase *>	vcfs { vcf_parents };
		new_merged_vcf = VCFGeno::join(vcfs, orig_vcf->get_samples());
	}
	delete vcf_parents;
	delete merged_vcf;
	return new_merged_vcf;
}
