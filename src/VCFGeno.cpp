#include "../include/VCFGeno.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFGenoBase ////////////////////

VCFGeno *VCFGenoBase::copy() const {
	vector<string>	new_samples(samples.size());
	std::copy(samples.begin(), samples.end(), new_samples.begin());
	vector<GenoRecord *>	new_records(size());
	for(size_t i = 0; i < size(); ++i) {
		new_records[i] = get_record(i)->copy();
	}
	return new VCFGeno(new_samples, new_records, vcf);
}

vector<GenoRecord *> VCFGenoBase::get_geno_records() const {
	vector<GenoRecord *>	records;
	for(size_t i = 0; i < size(); ++i)
		records.push_back(get_record(i));
	return records;
}

vector<size_t> VCFGenoBase::extract_columns(const STRVEC& samples) const {
	map<string, size_t>	dic;	// { sample: column }
	for(size_t i = 0; i < this->samples.size(); ++i)
		dic[this->samples[i]] = i;
	
	vector<size_t>	columns;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			columns.push_back(q->second);
		else
			columns.push_back(string::npos);
	}
	return columns;
}

vector<int> VCFGenoBase::extract_sample_genotypes(size_t sample_index) const {
	vector<int>	gts;
	for(size_t j = 0; j < size(); ++j) {
		const GenoRecord	*record = get_record(j);
		const int	gt = record->get_geno(sample_index);
		gts.push_back(gt);
	}
	return gts;
}

vector<int> VCFGenoBase::clip_raw_haplotype(size_t sample_index,
														int side) const {
	vector<int>	hap;
	for(size_t j = 0; j < size(); ++j) {
		const GenoRecord	*record = get_record(j);
		const int	gt = record->get_geno(sample_index);
		hap.push_back((gt >> side) & 1);
	}
	return hap;
}

vector<size_t> VCFGenoBase::map_samples_to_reference_columns() const {
	// Build a dictionary: reference sample name -> column index in VCF
	// (+9 because VCF sample columns start after the first 9 fixed fields)
	map<string, size_t>	dic_col;
	const STRVEC&	ref_samples = vcf->get_samples();
	for(size_t i = 0; i != ref_samples.size(); ++i) {
		dic_col[ref_samples[i]] = i + 9;
	}
	
	// Prepare result vector initialized with 0
	// (0 means the sample does not exist in the reference VCF)
	vector<size_t>	columns(samples.size(), 0);
	// For each sample in this VCF, look up its column index in the reference
	for(size_t i = 0; i < samples.size(); ++i) {
		auto	q = dic_col.find(samples[i]);
		if(q != dic_col.end())
			columns[i] = q->second;
	}
	return columns;
}

void VCFGenoBase::write(ostream& os, bool with_header) const {
	if(with_header) {
		const auto	header = vcf->trim_header(samples);
		for(auto p = header.begin(); p != header.end(); ++p)
			Common::write_tsv(*p, os);
	}
	
	const vector<size_t>	columns = map_samples_to_reference_columns();
	for(size_t i = 0; i < this->size(); ++i) {
		const GenoRecord	*record = this->get_record(i);
		record->write(vcf->get_record(i), columns, os);
	}
}

VCFGeno *VCFGenoBase::extract_by_samples(const STRVEC& samples) const {
	const auto	cs = extract_columns(samples);
	vector<GenoRecord *>	new_records;
	for(size_t i = 0; i < this->size(); ++i) {
		GenoRecord	*record = this->get_record(i);
		const ll	pos = record->get_pos();
		vector<int>	geno;
		for(auto p = cs.begin(); p != cs.end(); ++p) {
			geno.push_back(record->get_geno(*p));
		}
		GenoRecord	*new_record = new GenoRecord(pos, geno);
		new_records.push_back(new_record);
	}
	return new VCFGeno(samples, new_records, this->vcf);
}


//////////////////// VCFGeno ////////////////////

VCFGeno::~VCFGeno() {
	Common::delete_all(records);
}

vector<vector<int>> VCFGeno::create_ref_haps() const {
	const size_t	 N = num_samples() * 2;
	vector<vector<int>>	ref_haps(N, vector<int>(this->size(), 0));
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < this->size(); ++j)
			ref_haps[i][j] = (records[j]->get_geno(i>>1) >> (i&1)) & 1;
	}
	return ref_haps;
}

VCFGeno *VCFGeno::convert(const VCFSmall *vcf) {
	vector<GenoRecord *>	new_records;
	for(size_t i = 0; i < vcf->size(); ++i) {
		const VCFRecord	*record = vcf->get_record(i);
		const auto&	v = record->get_v();
		const ll	pos = record->pos();
		vector<int>	new_geno;
		for(auto p = v.begin()+9; p != v.end(); ++p)
			new_geno.push_back(Genotype::all_gt_to_int(*p));
		auto	*new_record = new GenoRecord(pos, new_geno);
		new_records.push_back(new_record);
	}
	return new VCFGeno(vcf->get_samples(), new_records, vcf);
}

VCFGeno *VCFGeno::join_core(const vector<const VCFGenoBase *>& vcfs,
							const SampleIndexMap& dic, const STRVEC& samples) {
	vector<tuple<string, const VCFGenoBase *, size_t>>	cols;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			cols.push_back(make_tuple(*p, q->second.first, q->second.second));
	}
	
	STRVEC	new_samples;
	for(auto p = cols.begin(); p != cols.end(); ++p)
		new_samples.push_back(get<0>(*p));
	
	vector<GenoRecord *>	new_records;
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		const ll	pos = vcfs[0]->get_record(i)->get_pos();
		vector<int>	geno;
		for(auto p = cols.begin(); p != cols.end(); ++p) {
			const auto	*vcf = get<1>(*p);
			const size_t	k = get<2>(*p);
			const auto	*record = vcf->get_record(i);
			const int	gt = record->get_geno(k);
			geno.push_back(gt);
		}
		new_records.push_back(new GenoRecord(pos, geno));
	}
	const VCFSmall	*ref_vcf = vcfs[0]->get_ref_vcf();
	return new VCFGeno(new_samples, new_records, ref_vcf);
}

VCFGeno *VCFGeno::join(const vector<const VCFGenoBase *>& vcfs,
											const STRVEC& samples) {
	// Register all samples with their corresponding VCF and sample index.
	SampleIndexMap	dic;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const VCFGenoBase	*vcf = *p;
		const STRVEC&	ss = vcf->get_samples();
		for(size_t i = 0; i < ss.size(); ++i)
			dic[ss[i]] = make_pair(vcf, i);
	}
	
	return VCFGeno::join_core(vcfs, dic, samples);
}

VCFGeno *VCFGeno::join(const VCFGenoBase *vcf1, const VCFGenoBase *vcf2,
														const STRVEC& samples) {
	if(vcf1 == NULL)
		return vcf2->copy();
	
	vector<const VCFGenoBase *>	vcfs { vcf1, vcf2 };
	return VCFGeno::join(vcfs, samples);
}

VCFGeno *VCFGeno::create_by_two_vcfs(const VCFGenoBase *vcf1,
											const VCFSmall *vcf2,
											const STRVEC& samples) {
	// assume that samples are [mat, pat, prog1, prog2, ...]
	const vector<size_t>	columns1 = vcf1->extract_columns(samples);
	const vector<size_t>	columns2 = vcf2->extract_columns(samples);
	vector<GenoRecord *>	new_records;
	for(size_t i = 0; i < vcf1->size(); ++i) {
		const auto	*record1 = vcf1->get_record(i);
		const auto	*record2 = vcf2->get_record(i);
		vector<int>	geno;
		for(size_t j = 0; j < samples.size(); ++j) {
			if(columns1[j] != string::npos)
				geno.push_back(record1->get_geno(columns1[j]));
			else if(columns2[j] != string::npos)
				geno.push_back(
					Genotype::all_gt_to_int(record2->get_gt(columns2[j]-9)));
			else
				geno.push_back(Genotype::NA);
		}
		auto	*new_record = new GenoRecord(record1->get_pos(), geno);
		new_records.push_back(new_record);
	}
	return new VCFGeno(samples, new_records, vcf2);
}

VCFGeno *VCFGeno::extract_samples(const STRVEC& samples, const VCFSmall *vcf) {
	const auto	cs = vcf->extract_columns(samples);
	vector<GenoRecord *>	new_records;
	for(size_t i = 0; i < vcf->size(); ++i) {
		const VCFRecord	*record = vcf->get_record(i);
		const ll	pos = record->pos();
		vector<int>	new_geno;
		for(auto p = cs.begin(); p != cs.end(); ++p) {
			const auto	c = *p;
			const int	gt = c == string::npos ? Genotype::NA :
									Genotype::all_gt_to_int(record->get_v()[c]);
			new_geno.push_back(gt);
		}
		auto	*new_record = new GenoRecord(pos, new_geno);
		new_records.push_back(new_record);
	}
	return new VCFGeno(samples, new_records, vcf);
}
