#include <sstream>
#include <cmath>
#include <numeric>
#include "../include/VCF.h"
#include "../include/common.h"

using namespace std;


//////////////////// VCFRecord ////////////////////

pair<string,ll> VCFRecord::position() const {
	return pair<string,ll>(v[0], stoll(v[1]));
}

string VCFRecord::gt(const string& sample) const {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p == sample)
			return v[p-samples.begin()+9];
	}
	return "";
}

STRVEC VCFRecord::gts() const {
	STRVEC	w(samples.size());
	std::copy(v.begin() + 9, v.end(), w.begin());
	return w;
}

#include <cassert>
vector<int> VCFRecord::get_int_gts() const {
	vector<int>	w(samples.size());
	for(size_t i = 0; i < samples.size(); ++i) {
		assert(i < w.size());
		w[i] = get_int_gt(i);
	}
	return w;
}

bool VCFRecord::is_homo(size_t i) const {
	const string&	s = this->v[i+9];
	const char	*ptr = s.c_str();
	return ptr[0] != '.' && ptr[0] == ptr[2];
}

bool VCFRecord::is_hetero(size_t i) const {
	const string&	s = this->v[i+9];
	const char	*ptr = s.c_str();
	return ptr[0] != '.' && ptr[2] != '.' && ptr[0] != ptr[2];
}

void VCFRecord::write(ostream& os) const {
	Common::write_tsv(this->v, os);
}

void VCFRecord::set_GT(size_t i, const string& gt) {
	this->v[i+9].replace(0, 3, gt);
}

void VCFRecord::set_int_GT(size_t i, int gt) {
	switch(gt) {
		case  0: set_GT(i, "0/0"); return;
		case  1: set_GT(i, "0/1"); return;
		case  2: set_GT(i, "1/1"); return;
		default: set_GT(i, "./."); return;
	}
}

WrongRecordType VCFRecord::check() const {
	if(v.size() != samples.size() + 9)
		return WrongRecordType::NUMCOLUMNERROR;
	
	for(auto p = v.begin() + 9; p != v.end(); ++p) {
		const string&	gt = *p;
		if(gt.size() < 3 || (gt.c_str()[1] != '/' && gt.c_str()[1] != '|'))
			return WrongRecordType::GENOTYPEERROR;
	}
	return WrongRecordType::RIGHT;
}

size_t VCFRecord::find_key_position(const string& key) const {
	return Genotype::find_key_position(v[8], key);
}

vector<VCFRecord::Probs> VCFRecord::parse_PL() const {
	vector<Probs>	probs;
	const size_t	N = this->num_samples();
	const size_t	PL_pos = this->find_key_position("PL");
	for(size_t i = 0; i < N; ++i) {
		const auto	w = Common::split(this->v[i+9], ':');
		try {
			// If PL_pos is std::string::npos,
			// throw std::out_of_range and
			// execute the processing in the catch block.
			const auto	pls = Common::split(w.at(PL_pos), ',');
			vector<double>	ps(3);
			for(size_t k = 0; k < 3; ++k) {
				ps[k] = exp(-log(10.0)*stoi(pls[k]));
			}
			const double	sum_ps = std::accumulate(ps.begin(), ps.end(), 0.0);
			probs.push_back(make_tuple(ps[0]/sum_ps,
										ps[1]/sum_ps, ps[2]/sum_ps));
		}
		catch(...) {
			switch(get_int_gt(i)) {
				case 0:  probs.push_back(make_tuple(1.0, 0.0, 0.0)); break;
				case 1:  probs.push_back(make_tuple(0.0, 1.0, 0.0)); break;
				case 2:  probs.push_back(make_tuple(0.0, 0.0, 1.0)); break;
				default: probs.push_back(make_tuple(1./3, 1./3, 1./3)); break;
			}
		}
	}
	return probs;
}


//////////////////// VCFBase ////////////////////

VCFBase::VCFBase(const vector<STRVEC>& h, const STRVEC& s) : 
								header(h), samples(s) {
	determine_chromosome_id();
}

POSITION VCFBase::position(const pair<string,ll>& p) const {
	const string&	chr = p.first;
	const ll		pos = p.second;
	auto	q = this->chrs.find(chr);
	if(q == this->chrs.end()) {
		this->chrs[chr] = (int)this->chrs.size() + 1;
	}
	return make_pair(this->chrs[chr], pos);
}

POSITION VCFBase::record_position(const VCFRecord& record) const {
	return this->position(record.position());
}

std::string VCFBase::chr(int chr_id) const {
	for(auto p = this->chrs.begin(); p != this->chrs.end(); ++p) {
		if(p->second == chr_id)
			return p->first;
	}
	return "";
}

int VCFBase::find_column(const string& sample) const {
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		if(*p == sample)
			return (p - samples.begin()) + 9;
	}
	return -1;
}

void VCFBase::write_header(ostream& os) const {
	for(auto p = header.begin(); p != header.end(); ++p)
		Common::write_tsv(*p, os);
}

void VCFBase::determine_chromosome_id() {
	int	id = 0;
	for(auto p = header.begin(); p != header.end(); ++p) {
		const string&	s = p->front();
		if(s.find("##config=<ID=") == string::npos)
			continue;
		size_t	q = s.find(",length");
		string	chr	= s.substr(13U, q);
		chrs[chr] = id;
		id += 1;
	}
}


//////////////////// VCFReader ////////////////////

VCFReader::VCFReader(const string& path) {
	if(path.substr(path.length() - 3U) != ".gz")
		reader = new FileReader(path);
	else
		reader = new FileReaderGZ(path);
}

VCFReader::~VCFReader() {
	delete reader;
}

void VCFReader::read_header() {
	string	line;
	while(reader->getline(line)) {
		header.push_back(Common::split(line, '\t'));
		if(line.substr(0U, 6U) == "#CHROM")
			break;
	}
}

STRVEC VCFReader::next() {
	string	line;
	if(!reader->getline(line))
		return STRVEC();
	
	return Common::split(line, '\t');
}

STRVEC VCFReader::get_samples() const {
	const STRVEC&	v = header.back();
	STRVEC	samples(v.begin() + 9, v.end());
	return samples;
}


//////////////////// VCFSmallBase ////////////////////

vector<STRVEC> VCFSmallBase::trim_header(const STRVEC& samples) const {
	vector<STRVEC>	header = this->get_header();
	STRVEC	v(header.back().begin(), header.back().begin() + 9);
	v.insert(v.end(), samples.begin(), samples.end());
	header.back() = v;
	return header;
}

vector<int> VCFSmallBase::clip_raw_haplotype(size_t sample_id, int i) const {
	vector<int>	hap;
	for(size_t j = 0; j < size(); ++j) {
		const VCFRecord	*record = get_record(j);
		const string&	gt = record->get_gt(sample_id);
		const char	c = gt.c_str()[i*2];
		hap.push_back(c == '.' ? -1 : (int)(c - '0'));
	}
	return hap;
}

vector<size_t> VCFSmallBase::extract_columns(STRVEC::const_iterator first,
											STRVEC::const_iterator last) const {
	const STRVEC&	samples = this->get_samples();
	map<string, size_t>	dic;
	for(size_t i = 0; i < samples.size(); ++i)
		dic[samples[i]] = i + 9;
	
	vector<size_t>	columns;
	for(auto p = first; p != last; ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			columns.push_back(q->second);
		else
			columns.push_back(string::npos);
	}
	return columns;
}

VCFSmall *VCFSmallBase::extract_samples(const STRVEC& samples) const {
	const auto	header = trim_header(samples);
	vector<VCFRecord *>	empty_records;
	VCFSmall	*vcf = new VCFSmall(header, samples, empty_records);
	const vector<size_t>	cs = extract_columns(samples);
	for(size_t i = 0; i < size(); ++i) {
		const VCFRecord	*record = get_record(i);
		const STRVEC&	v = record->get_v();
		STRVEC	new_v(v.begin(), v.begin() + 9);
		for(auto q = cs.begin(); q != cs.end(); ++q)
			new_v.push_back(v[*q]);
		VCFRecord	*new_record = new VCFRecord(new_v, vcf->get_samples());
		vcf->add_record(new_record);
	}
	return vcf;
}

void VCFSmallBase::write(ostream& os, bool write_header) const {
	if(write_header)
		this->write_header(os);
	for(size_t i = 0; i < size(); ++i)
		get_record(i)->write(os);
}

void VCFSmallBase::write_header(ostream& os) const {
	const auto&	header = get_header();
	for(auto p = header.begin(); p != header.end(); ++p)
		Common::write_tsv(*p, os);
}


//////////////////// VCFSmall ////////////////////

VCFSmall::VCFSmall(const vector<STRVEC>& h, const STRVEC& s,
								const vector<VCFRecord *>& rs, bool rr) :
				VCFBase(h, s), VCFSmallBase(), records(rs), reuses_records(rr) {
	for(auto p = records.begin(); p != records.end(); ++p)
		this->record_position(**p);
}

VCFSmall::~VCFSmall() {
	if(!reuses_records) {
		for(auto p = records.begin(); p != records.end(); ++p)
			delete *p;
	}
}

void VCFSmall::check_records() const {
	int	wrong_counter = 0;
	vector<pair<VCFRecord *, WrongRecordType>>	wrong_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const WrongRecordType	type = (*p)->check();
		if(type != WrongRecordType::RIGHT) {
			wrong_counter += 1;
			// Display errors for up to five records
			if(wrong_counter <= 5)
				wrong_records.push_back(make_pair(*p, type));
		}
	}
	if(wrong_counter > 0) {
		throw RecordException(wrong_counter, wrong_records);
	}
}

VCFSmall *VCFSmall::select_samples(const vector<string>& new_samples) const {
	const size_t	N = new_samples.size();
	map<string, size_t>	dic;
	for(size_t i = 0; i != samples.size(); ++i) {
		dic[samples[i]] = i + 9;
	}
	
	vector<size_t>	cs;
	for(auto p = new_samples.begin(); p != new_samples.end(); ++p) {
		cs.push_back(dic[*p]);
	}
	
	vector<VCFRecord *>	new_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		vector<string>	v(N+9);
		std::copy((*p)->get_v().begin(), (*p)->get_v().begin()+9, v.begin());
		for(size_t i = 0; i < N; ++i) {
			v[i+9] = (*p)->get_v()[cs[i]];
		}
		new_records.push_back(new VCFRecord(v, new_samples));
	}
	const auto	header = trim_header(new_samples);
	return new VCFSmall(header, new_samples, new_records);
}

VCFSmall *VCFSmall::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	
	vector<VCFRecord *>	records;
	while(true) {
		const STRVEC	v = reader->next();
		if(v.empty())
			break;
		VCFRecord	*record = new VCFRecord(v, samples);
		records.push_back(record);
	}
	delete reader;
	
	return new VCFSmall(header, samples, records);
}

// join VCFs in order of given samples
VCFSmall *VCFSmall::join(const vector<const VCFSmallBase *>& vcfs,
													const STRVEC& samples) {
	map<string, pair<const VCFSmallBase *, size_t>>	dic;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const VCFSmallBase	*vcf = *p;
		const STRVEC&	ss = vcf->get_samples();
		for(size_t i = 0; i < ss.size(); ++i)
			dic[ss[i]] = make_pair(vcf, i + 9);
	}
	
	vector<tuple<string, const VCFSmallBase *, size_t>>	cols;
	for(auto p = samples.begin(); p != samples.end(); ++p) {
		auto	q = dic.find(*p);
		if(q != dic.end())
			cols.push_back(make_tuple(*p, q->second.first, q->second.second));
	}
	
	STRVEC	new_samples;
	for(auto p = cols.begin(); p != cols.end(); ++p)
		new_samples.push_back(get<0>(*p));
	
	// new_vcf owns samples
	auto	new_header = vcfs.front()->trim_header(new_samples);
	vector<VCFRecord *>	empty_records;
	VCFSmall	*new_vcf = new VCFSmall(new_header, new_samples, empty_records);
	const STRVEC&	samples_ = new_vcf->get_samples();
	for(size_t i = 0; i < vcfs.front()->size(); ++i) {
		const VCFSmallBase	*vcf = vcfs[0];
		const STRVEC&	orig_v = vcf->get_record(i)->get_v();
		STRVEC	v(orig_v.begin(), orig_v.begin() + 9);
		for(auto p = cols.begin(); p != cols.end(); ++p)
			v.push_back(get<1>(*p)->get_record(i)->get_v()[get<2>(*p)]);
		new_vcf->add_record(new VCFRecord(v, samples_));
	}
	
	return new_vcf;
}

VCFSmall *VCFSmall::join(const VCFSmallBase *vcf1, const VCFSmallBase *vcf2,
														const STRVEC& samples) {
	vector<const VCFSmallBase *>	vcfs{ vcf1, vcf2 };
	return VCFSmall::join(vcfs, samples);
}

VCFSmall *VCFSmall::create_by_two_vcfs(const VCFSmallBase *vcf1,
										 const VCFSmallBase *vcf2,
										 const STRVEC& samples) {
	// assume that samples are [mat, pat, prog1, prog2, ...]
	const vector<size_t>	columns1 = vcf1->extract_columns(samples);
	const vector<size_t>	columns2 = vcf2->extract_columns(samples);
	const auto	new_header = vcf1->trim_header(samples);
	
	vector<VCFRecord *>	new_records;
	for(size_t i = 0; i < vcf1->size(); ++i) {
		const auto	*record1 = vcf1->get_record(i);
		const auto	*record2 = vcf2->get_record(i);
		STRVEC	v(record1->get_v().begin(), record1->get_v().begin() + 9);
		for(size_t j = 0; j < samples.size(); ++j) {
			if(columns1[j] != string::npos)
				v.push_back(record1->get_v()[columns1[j]]);
			else if(columns2[j] != string::npos)
				v.push_back(record2->get_v()[columns2[j]]);
			else
				v.push_back("./.");
		}
		auto	*new_record = new VCFRecord(v, samples);
		new_records.push_back(new_record);
	}
	return new VCFSmall(new_header, samples, new_records);
}


//////////////////// VCFHuge::ChromDivisor ////////////////////

VCFHuge::ChromDivisor::~ChromDivisor() {
	Common::delete_all(this->records);
}

VCFSmall *VCFHuge::ChromDivisor::next() {
	while(true) {
		VCFRecord	*record = this->vcf->next();
		if(state == STATE::START) {
			if(record == NULL)
				return NULL;
			else {
				this->chrom = record->chrom();
				this->records.push_back(record);
				state = STATE::DOING;
			}
		}
		else if(state == STATE::END) {
			return NULL;
		}
		else if(record == NULL) {
			VCFSmall	*vcf_chrom = new VCFSmall(this->vcf->header,
											this->vcf->samples, this->records);
			this->records.clear();
			state = STATE::END;
			return vcf_chrom;
		}
		else if(record->chrom() != this->chrom) {
			VCFSmall	*vcf_chrom = new VCFSmall(this->vcf->header,
											this->vcf->samples, this->records);
			this->records.clear();
			this->records.push_back(record);
			this->chrom = record->chrom();
			return vcf_chrom;
		}
		else {
			this->records.push_back(record);
		}
	}
}


//////////////////// VCFHuge ////////////////////

VCFHuge::VCFHuge(const vector<STRVEC>& h, const STRVEC& s, VCFReader *r_) :
												VCFBase(h, s), reader(r_) { }

VCFHuge::~VCFHuge() {
	delete this->reader;
}

void VCFHuge::leave_only_GT(vector<string>& v) {
	v[8] = "GT";
	for(auto p = v.begin() + 9; p != v.end(); ++p) {
		const size_t	q = p->find(':');
		*p = p->substr(0, q);
	}
}

VCFRecord *VCFHuge::next() {
	STRVEC	v = reader->next();
	if(v.empty())
		return NULL;
	
//	leave_only_GT(v);
	VCFRecord	*record = new VCFRecord(v, this->samples);
	this->record_position(*record);		// necessary to update chrs
	return record;
}

VCFRecord *VCFHuge::proceed(POSITION pos) {
	while(true) {
		VCFRecord	*record = next();
		if(record == NULL)
			return NULL;
		
		if(this->record_position(*record) == pos)
			return record;
		else
			delete record;
	}
	return NULL;
}

VCFHuge *VCFHuge::read(const string& path) {
	VCFReader	*reader = new VCFReader(path);
	reader->read_header();
	const vector<STRVEC>&	header = reader->get_header();
	const STRVEC	samples = reader->get_samples();
	return new VCFHuge(header, samples, reader);
}


//////////////////// RecordException ////////////////////

RecordException::RecordException(int counter,
						const vector<pair<VCFRecord *, WrongRecordType>>& rs) {
	stringstream	ss;
	if(counter == 1)
		ss << "error : 1 record is wrong :";
	else
		ss << "error : " << counter << " records are wrong :";
	
	for(auto p = rs.begin(); p != rs.end(); ++p) {
		const auto&	v = p->first->get_v();
		ss << '\n' << v[0] << '\t' << v[1] << " : ";
		if(p->second == WrongRecordType::GENOTYPEERROR)
			ss << "wrong genotype.";
		else
			ss << "the number of items is wrong.";
	}
	
	message = ss.str();
}

const char *RecordException::what() const noexcept {
	return message.c_str();
}
