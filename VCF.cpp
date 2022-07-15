#include <algorithm>
#include "VCF.h"
#include "common.h"

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

int VCFRecord::get_int_gt(size_t i) const {
	const string&	s = v[i+9];
	try {
		return stoi(s.substr(0U, 1U)) + stoi(s.substr(2U, 1U));
	}
	catch(std::invalid_argument& e) {
		return -1;
	}
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

STRVEC VCFRecord::extract_v(const STRVEC& samples) const {
	map<string,size_t>	dic;
	for(size_t i = 0U; i < this->samples.size(); ++i)
		dic[this->samples[i]] = i;
	
	STRVEC	new_v(v.begin(), v.begin() + 9);
	for(auto p = samples.begin(); p != samples.end(); ++p)
		new_v.push_back(v[dic[*p]+9]);
	return new_v;
}

void VCFRecord::write(ostream& os) const {
	Common::write_tsv(this->v, os);
}

void VCFRecord::copy_properties(STRVEC::iterator it) const {
	std::copy(v.begin(), v.begin() + 9, it);
}

void VCFRecord::set_GT(size_t i, const string& gt) {
	this->v[i].replace(0, 3, gt);
}

void VCFRecord::set_GTs(const STRVEC& GTs) {
	for(size_t i = 9U; i < v.size(); ++i)
		this->set_GT(i, GTs[i]);
}

void VCFRecord::set_int_GT(size_t i, int gt) {
	switch(gt) {
		case  0: set_GT(i, "0/0"); return;
		case  1: set_GT(i, "0/1"); return;
		case  2: set_GT(i, "1/1"); return;
		default: set_GT(i, "./."); return;
	}
}


//////////////////// VCFBase ////////////////////

VCFBase::VCFBase(const vector<STRVEC>& h, const STRVEC& s) : 
								header(h), samples(s),
								sample_ids(number_samples(s)) {
	determine_chromosome_id();
}

map<string,size_t> VCFBase::number_samples(const STRVEC& samples_) const {
	map<string,size_t>	dic;
	size_t	id = 0U;
	for(auto p = samples_.begin(); p != samples_.end(); ++p) {
		dic[*p] = id;
		++id;
	}
	return dic;
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


//////////////////// VCFHuge ////////////////////

VCFHuge::VCFHuge(const vector<STRVEC>& h, const STRVEC& s, VCFReader *r_) :
												VCFBase(h, s), reader(r_) { }

VCFHuge::~VCFHuge() {
	delete reader;
}

VCFRecord *VCFHuge::next() {
	const STRVEC	v = reader->next();
	if(v.empty())
		return NULL;
	
	VCFRecord	*record = new VCFRecord(v, this->samples);
	this->record_position(*record);		// chrsを更新するために必要
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


//////////////////// VCFSmall ////////////////////

VCFSmall::VCFSmall(const vector<STRVEC>& h, const STRVEC& s,
											vector<VCFRecord *> rs) :
												VCFBase(h, s), records(rs) {
	for(auto p = records.begin(); p != records.end(); ++p)
		this->record_position(**p);
}

VCFSmall::~VCFSmall() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

void VCFSmall::write(ostream& os) const {
	this->write_header(os);
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->write(os);
}

void VCFSmall::update_genotypes(const std::vector<STRVEC>& GT_table) {
	for(size_t i = 0U; i < records.size(); ++i)
		records[i]->set_GTs(GT_table[i]);
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
