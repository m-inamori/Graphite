#include <sstream>
#include <stack>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/VCFSelfHetero.h"
#include "../include/common.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"
#include "../include/inverse_graph.h"
#include "../include/option.h"

using namespace std;


//////////////////// VCFSelfHeteroRecord ////////////////////

vector<int> VCFSelfHeteroRecord::progeny_gts() const {
	const vector<int>	gts = this->get_int_gts();
	return vector<int>(gts.begin() + 1, gts.end());
}

void VCFSelfHeteroRecord::set_haplo(int h) {
	this->set_GT(0, h == 0 ? "0|1" : "1|0");
}

void VCFSelfHeteroRecord::set_int_gt_by_which_comes_from(int w1, int w2,
																	size_t i) {
	const string&	parent_gt = this->get_gt(0);
	stringstream	ss;
	ss << parent_gt.c_str()[w1*2] << "|" << parent_gt.c_str()[w2*2];
	const string&	GT = ss.str();
	set_GT(i + 1, GT);
}


//////////////////// VCFSelfHetero ////////////////////

VCFSelfHetero::VCFSelfHetero(const vector<STRVEC>& h, const STRVEC& s,
								const vector<VCFSelfHeteroRecord *>& rs,
								const Map& m) :
							VCFBase(h, s), VCFSmallBase(),
							VCFMeasurable(m), records(rs) { }

VCFSelfHetero::~VCFSelfHetero() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

InvGraph VCFSelfHetero::make_graph(double max_dist) const {
	vector<vector<int>>	gtss;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		gtss.push_back((*p)->progeny_gts());
	}
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	// Graph creation is time-consuming with all record pairs.
	// Only pairs within 10cM are checked.
	// If records are scarce, the graph may not connect within this range.
	// All pairs are checked for up to 30 records.
	// For larger record counts,
	// the calculation is adjusted to be equivalent without considering cM.
	// However, at least 10 preceding and following records are always checked.
	const size_t	range = std::max(10UL, std::min(L, 900/L));
	
	InvGraph	graph;
	// isolated nodes can not be represented without pre-registering the key
	for(size_t k = 0U; k < L; ++k)
		graph[k];
	for(size_t k = 0; k < L; ++k) {
		for(size_t l = k + 1; l < L; ++l) {
			const double	cM = cMs[l] - cMs[k];
			// don't see a connection
			// if the records are more than 10 cM apart
			if(cM > 10.0 && k + range < l)
				break;
			const auto	p = distance(gtss[k], gtss[l]);
			const double	dist = p.first;
			const bool		inv = p.second;
			if(dist <= max_dist) {
				graph[k].push_back(make_tuple(l, dist, inv));
				graph[l].push_back(make_tuple(k, dist, inv));
			}
		}
	}
	return graph;
}

// Since beta distribution is used, dist is a double
pair<double, bool> VCFSelfHetero::distance(const vector<int>& gts1,
											const vector<int>& gts2) {
	const int	N = (int)gts1.size();
	int	counter_homo_right = 0;
	int	counter_hetero_right = 0;
	int	counter_inverse = 0;
	int	counter_NA = 0;
	for(size_t i = 0; i < gts1.size(); ++i) {
		const int	int_gt1 = gts1[i];
		const int	int_gt2 = gts2[i];
		if(int_gt1 == -1 || int_gt2 == -1)
			counter_NA += 1;
		else if(int_gt1 == 1 && int_gt2 == 1)
			counter_hetero_right += 1;
		else if(int_gt1 == int_gt2)
			counter_homo_right += 1;
		else if((int_gt1 == 0 && int_gt2 == 2) ||
				(int_gt1 == 2 && int_gt2 == 0))
			counter_inverse += 1;
	}
	
	if(counter_homo_right >= counter_inverse) {
		const int	right = counter_hetero_right + counter_homo_right;
		return make_pair(dist_with_NA(right, counter_NA, N), false);
	}
	else {
		const int	right = counter_hetero_right + counter_inverse;
		return make_pair(dist_with_NA(right, counter_NA, N), true);
	}
}

// beta distribution
double VCFSelfHetero::dist_with_NA(int right, int counter_NA, int N) {
	const int	diff = N - right - counter_NA;
	const double	diff_ratio = (diff + 1.0) / (right + diff + 2.0);
	return diff + counter_NA * diff_ratio;
}

// decide heterozygous parent haplotypes
vector<int> VCFSelfHetero::make_parent_haplotypes(const InvGraph& graph) const {
	const InvGraph	tree = graph.minimum_spanning_tree();
	auto	vs = tree.collect_nodes();
	std::sort(vs.begin(), vs.end());
	vector<int>	haplo = create_haplotype(vs[0], tree);
	
	vector<int>	haplo2;
	for(auto p = haplo.begin(); p != haplo.end(); ++p) {
		if(*p != -1)
				haplo2.push_back(*p);
	}
	return haplo2;
}

vector<int> VCFSelfHetero::create_haplotype(size_t v0,
											const InvGraph& tree) const {
	const size_t	L = size();
	vector<int>	haplo(L, -1);
	haplo[v0] = 0;
	const auto	edges = tree.walk(v0);
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const size_t	v1 = get<0>(*p);
		const size_t	v2 = get<1>(*p);
		const bool		inv = get<2>(*p);
		if(inv)
			haplo[v2] = haplo[v1] == 0 ? 1 : 0;		// inverse
		else
			haplo[v2] = haplo[v1];
	}
	return haplo;
}

VCFSelfHetero *VCFSelfHetero::make_subvcf(const InvGraph& graph) const {
	const vector<int>	haplo = make_parent_haplotypes(graph);
	vector<int>	indices;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		indices.push_back(p->first);
	std::sort(indices.begin(), indices.end());
	vector<VCFSelfHeteroRecord *>	records;
	for(size_t i = 0; i < haplo.size(); ++i) {
		VCFSelfHeteroRecord	*record = this->records[indices[i]];
		record->set_haplo(haplo[i]);
		records.push_back(record);
	}
	return new VCFSelfHetero(this->header, this->samples,
										records, this->get_map());
}

pair<vector<VCFSelfHetero *>, vector<VCFSelfHeteroRecord *>>
		VCFSelfHetero::determine_haplotype(const OptionImpute *option) const {
	const double	max_dist = std::min((double)option->max_dist,
										this->num_progenies() * 0.1) * 2;
	const InvGraph	graph = this->make_graph(max_dist);
	const vector<InvGraph>	subgraphs = graph.divide_into_connected();
	
	// keep only large graphs and collect records of small graphs
	vector<InvGraph>	large_graphs;
	vector<InvGraph>	small_graphs;
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		if((int)p->size() >= option->min_graph)
			large_graphs.push_back(*p);
		else
			small_graphs.push_back(*p);
	}
	
	// keep only the largest graph if there are only small ones
	if(large_graphs.empty()) {
		size_t	max_i = 0;
		for(size_t i = 1; i < small_graphs.size(); ++i) {
			if(small_graphs[i].size() > small_graphs[max_i].size())
				max_i = i;
		}
		large_graphs.push_back(small_graphs[max_i]);
		small_graphs.erase(small_graphs.begin() + max_i);
	}
	
	vector<VCFSelfHetero *>	subvcfs;
	for(auto p = large_graphs.begin(); p != large_graphs.end(); ++p) {
		subvcfs.push_back(this->make_subvcf(*p));
	}
	
	vector<VCFSelfHeteroRecord *>	unused_records;
	for(auto p = small_graphs.begin(); p != small_graphs.end(); ++p) {
		for(auto q = p->begin(); q != p->end(); ++q)
			unused_records.push_back(this->records[q->first]);
	}
	
	return make_pair(subvcfs, unused_records);
}

string VCFSelfHetero::make_seq(size_t i) const {
	string	seq;
	const size_t	q = i >> 1;
	const size_t	r = i & 1;
	for(size_t j = 0; j < this->size(); ++j) {
		const VCFSelfHeteroRecord	*record = records[j];
		const char	a = record->get_gt(q+1).c_str()[r*2];
		if(a == record->get_gt(0).c_str()[0])
			seq.push_back('0');
		else if(a != '.')
			seq.push_back('1');
		else
			seq.push_back('N');
	}
	return seq;
}

string VCFSelfHetero::impute_each_sample_seq(size_t i,
								const vector<double>& cMs, double min_c) {
	const string	seq = this->make_seq(i);
	if(Imputer::is_all_same_without_N(seq))
		return Imputer::create_same_color_string(seq, '0');
	
	const string	hidden_seq = Imputer::impute(seq, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

void VCFSelfHetero::impute_each_sample(size_t i, const vector<double>& cMs,
																double min_c) {
	const string	imputed_seq1 = impute_each_sample_seq(i*2, cMs, min_c);
	const string	imputed_seq2 = impute_each_sample_seq(i*2+1, cMs, min_c);
	for(size_t k = 0; k < this->size(); ++k) {
		VCFSelfHeteroRecord	*record = this->records[k];
		const int	w1 = imputed_seq1.c_str()[k] - '0';
		const int	w2 = imputed_seq2.c_str()[k] - '0';
		record->set_int_gt_by_which_comes_from(w1, w2, i);
	}
}

void VCFSelfHetero::impute_each(const OptionImpute *option) {
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t k = 0; k < L; ++k) {
		cMs.push_back(this->record_cM(k));
	}
	
	const double	min_c = option->min_crossover;
	for(size_t i = 0; i < num_progenies(); ++i) {
		impute_each_sample(i, cMs, min_c);
	}
}

const OptionImpute *VCFSelfHetero::create_option(int num_threads) const {
	const size_t	num = std::max(2UL, this->size());
	const double	cM_length = this->record_cM(records.size()-1);
	const int	max_dist = std::max(4,
						(int)(cM_length * num_progenies()
										/ num * log10(num) * 2.5 * 0.01));
	return new OptionImpute(max_dist, 20, 5, 1.0, num_threads);
}

pair<vector<VCFSelfHetero *>, vector<VCFSelfHeteroRecord *>>
										VCFSelfHetero::impute(int num_threads) {
	if(this->records.empty()) {
		// If it do not create a new VCF,
		// it will delete the VCF that cannot be deleted.
		VCFSelfHetero	*empty_vcf = new VCFSelfHetero(header, samples,
								vector<VCFSelfHeteroRecord *>(), get_map());
		return make_pair(vector<VCFSelfHetero *>(1, empty_vcf),
									vector<VCFSelfHeteroRecord *>());
	}
	
	const OptionImpute	*option = this->create_option(num_threads);
	auto	p = this->determine_haplotype(option);
	vector<VCFSelfHetero *>&	vcfs = p.first;
	vector<VCFSelfHeteroRecord *>&	unused_records = p.second;
	std::for_each(unused_records.begin(), unused_records.end(),
					std::mem_fun(&VCFImpSelfRecord::enable_modification));
	for(auto q = vcfs.begin(); q != vcfs.end(); ++q)
		(*q)->impute_each(option);
	delete option;
	return make_pair(vcfs, unused_records);
}
