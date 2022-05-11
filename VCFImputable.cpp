#include <sstream>
#include <set>
#include <algorithm>
#include <cassert>
#include <pthread.h>
#include "VCFImputable.h"
#include "Map.h"
#include "VCFCollection.h"
#include "paint.h"
#include "common.h"
#include "option.h"

using namespace std;


//////////////////// VCFImputableRecord ////////////////////

VCFImputableRecord::VCFImputableRecord(const STRVEC& v, const STRVEC& s,
										int h1, int h2, bool b, size_t i,
										const std::vector<int>& igs) :
									VCFHeteroHomoRecord(v, s), mat_hetero(b),
									int_gts(igs.empty() ? get_int_gts() : igs),
									index(i), group_id(0) {
	haplo[0] = h1;
	haplo[1] = h2;
}

STRVEC VCFImputableRecord::get_GTs() const {
	STRVEC	GTs;
	for(auto p = v.begin() + 9; p != v.end(); ++p)
		GTs.push_back(p->substr(0, 3));
	return GTs;
}

int VCFImputableRecord::difference(const VCFImputableRecord *other,
														bool inv) const {
	int diff = 0;
	for(size_t i = 0U; i < which_comes_from.size(); ++i) {
		if((which_comes_from[i] == other->which_comes_from[i]) == inv)
			diff += 1;
	}
	return diff;
}

int VCFImputableRecord::which_comes_from_hetero_parent_chromosomes_each(
																	int gt) {
	if(gt == -1)
		return -1;
	
	const int	gt_ = homo_int_gt() == 0 ? gt : gt - 1;
	if(gt_ != 0 && gt_ != 1)
		return -1;
	else if(haplo[0] == 0)
		return gt_;
	else
		return 1 - gt_;
}

void VCFImputableRecord::which_comes_from_hetero_parent_chromosomes() {
	for(size_t i = 2U; i < int_gts.size(); ++i) {
		const int	gt = int_gts[i];
		const int	w = which_comes_from_hetero_parent_chromosomes_each(gt);
		this->which_comes_from.push_back(w);
	}
}

int VCFImputableRecord::recover_int_gt(int h) {
	if(h == -1)
		return -1;
	else
		return haplo[h] + (homo_int_gt() >> 1);
}

void VCFImputableRecord::set_int_gts_by_which_comes_from(
											const vector<int>& hs) {
	this->which_comes_from = hs;
	for(size_t i = 2U; i < samples.size(); ++i) {
		int_gts[i] = recover_int_gt(hs[i-2]);
	}
}

int VCFImputableRecord::GT_position() const {
	const auto	w = Common::split(format(), ':');
	for(int i = 0; i < (int)w.size(); ++i) {
		if(w[i] == "GT")
			return i;
	}
	return -1;
}

string VCFImputableRecord::decode_gt(int h) {
	if(h == -1)
		return "./.";
	
	stringstream	ss;
	if(is_mat_hetero())
		ss << haplo[h] << '|' << (pat_int_gt() >> 1);
	else
		ss << (mat_int_gt() >> 1) << '|' << haplo[h];
	return ss.str();
}

string VCFImputableRecord::decode_parent_gt(int int_gt) {
	stringstream	ss;
	if(int_gt == 1)
		ss << haplo[0] << '|' << haplo[1];
	else
		ss << (int_gt >> 1) << '|' << (int_gt >> 1);
	return ss.str();
}

// ex) '1/1:10', '1|1' -> '1|1:10'
string VCFImputableRecord::update_gt(const string& gt,
										const string& s, int gt_pos) {
	auto	w = Common::split(s, ':');
	w[gt_pos] = gt;
	return Common::join(w, ':');
}

void VCFImputableRecord::update_genotypes() {
	const int	gt_pos = GT_position();
	assert(gt_pos != -1);
	v[9] = update_gt(decode_parent_gt(int_gts[0]), v[9], gt_pos);
	v[10] = update_gt(decode_parent_gt(int_gts[1]), v[10], gt_pos);
	for(size_t i = 11U; i < samples.size() + 9; ++i)
		v[i] = update_gt(decode_gt(which_comes_from[i-11]), v[i], gt_pos);
}

void VCFImputableRecord::inverse_haplotype() {
	const int	tmp = haplo[0];
	haplo[0] = haplo[1];
	haplo[1] = tmp;
	for(size_t i = 0U; i < which_comes_from.size(); ++i)
		which_comes_from[i] = which_comes_from[i] == 0 ? 1 : 0;
}


//////////////////// VCFImputable ////////////////////

VCFImputable::VCFImputable(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFImputableRecord *> rs, const Map& m) :
						VCFHeteroHomo(h, s, to_VCFHeteroHomoRecord(rs), m),
						records(rs),
						mat_hetero(determine_which_parents_is_hetero(this)) { }

vector<VCFHeteroHomoRecord *> VCFImputable::to_VCFHeteroHomoRecord(
											vector<VCFImputableRecord *>& rs) {
	vector<VCFHeteroHomoRecord *>	records(rs.size());
	std::copy(rs.begin(), rs.end(), records.begin());
	return records;
}

vector<STRVEC> VCFImputable::get_GT_table() const {
	vector<STRVEC>	table;
	for(auto p = records.begin(); p != records.end(); ++p)
		table.push_back((*p)->get_GTs());
	return table;
}

void VCFImputable::set_group_id(int id) {
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->set_group_id(id);
}

Graph::WeightedGraph VCFImputable::trim_inverse(
								const Graph::InvGraph& graph) const {
	Graph::WeightedGraph	weighted_graph;
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const auto&	v = p->second;
		for(auto q = v.begin(); q != v.end(); ++q)
			weighted_graph[p->first].push_back(
								pair<size_t,int>(get<0>(*q), get<1>(*q)));
	}
	return weighted_graph;
}

Graph::InvGraph VCFImputable::filter_graph(const Graph::InvGraph& graph,
									const Graph::WeightedGraph& tree) const {
	set<pair<size_t,size_t>>	edges;
	for(auto p = tree.begin(); p != tree.end(); ++p) {
		const size_t	v1 = p->first;
		const auto&	vs = p->second;
		for(auto q = vs.begin(); q != vs.end(); ++q)
			edges.insert(pair<size_t,size_t>(v1, q->first));
	}
	
	Graph::InvGraph	tree_graph;
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const size_t	v1 = p->first;
		const auto&		vs = p->second;
		for(auto q = vs.begin(); q != vs.end(); ++q) {
			const size_t	v2 = get<0>(*q);
			if(edges.find(pair<size_t,size_t>(v1, v2)) != edges.end())
				tree_graph[v1].push_back(*q);
		}
	}
	return tree_graph;
}

bool VCFImputable::is_all_same_without_N(const string& seq) const {
	char	c = '.';
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			if(c == '.')	// initial
				c = *p;
			else if(*p != c)
				return false;
		}
	}
	return true;
}

string VCFImputable::create_same_color_string(const string& seq) const {
	char	c = '0';	// dummy
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			c = *p;
			break;
		}
	}
	return std::string(seq.size(), c);
}

string VCFImputable::impute_each_seq(const string& seq,
										const vector<Matrix>& Ts,
										double MIN_CROSSOVER) const {
	if(is_all_same_without_N(seq))
		return create_same_color_string(seq);
	
	const vector<char>	hidden_states = { '0', '1' };
	const vector<char>	states = { '0', '1', 'N' };
	const string	hidden_seq = BaumWelch::impute(seq,
												hidden_states, states, Ts);
	
	const auto&	rs = VCFFamily::records;
	const size_t	L = rs.size();
	vector<double>	cMs(L);
	for(size_t k = 0U; k < L; ++k)
		cMs[k] = this->cM(k);
	const string	painted_seq = Paint::paint(hidden_seq, cMs, MIN_CROSSOVER);
	return painted_seq;
}

vector<VCFImputable::Matrix>
VCFImputable::create_transition_probability_matrix() const {
	const size_t	L = VCFFamily::size();
	vector<double>	cMs(L);
	for(size_t k = 0U; k < L; ++k)
		cMs[k] = this->cM(k);
	
	vector<Matrix>	Ts;
	for(auto q = cMs.begin(); q != cMs.end() - 1; ++q) {
		const double	p = (*(q+1) - *q) / 100;
		Matrix	T;
		T[pair<char,char>('0', '0')] = 1.0 - p;
		T[pair<char,char>('0', '1')] = p;
		T[pair<char,char>('1', '0')] = p;
		T[pair<char,char>('1', '1')] = 1.0 - p;
		Ts.push_back(T);
	}
	return Ts;
}

string VCFImputable::make_seq(size_t i) const {
	string	seq;
	for(auto p = records.begin(); p != records.end(); ++p) {
		const auto	*record = *p;
		const int	gt = record->get_which_comes_from(i);
		if(gt == 0)
			seq.push_back('0');
		else if(gt == 1)
			seq.push_back('1');
		else
			seq.push_back('N');
	}
	return seq;
}

void VCFImputable::impute(double MIN_CROSSOVER, int T) {
	vector<string>	hidden_seqs(samples.size() - 2);
	const auto	Ts = create_transition_probability_matrix();
	
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(this, Ts, i,
										MIN_CROSSOVER, T, hidden_seqs);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&impute_by_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		thread_function(configs[i]);
#endif
	
	for(int i = 0; i < T; ++i)
		delete configs[i];
	
	for(size_t i = 0U; i < records.size(); ++i) {
		auto	*record = records[i];
		vector<int>	hs;
		for(size_t j = 0U; j < VCFFamily::samples.size() - 2; ++j) {
			const int	h = hidden_seqs[j].c_str()[i] - '0';
			hs.push_back(h);
		}
		record->set_int_gts_by_which_comes_from(hs);
		record->update_genotypes();
	}
}

void VCFImputable::update_genotypes() {
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->update_genotypes();
	
	const auto	GT_table = get_GT_table();
	VCFHeteroHomo::update_genotypes(GT_table);	// 親も更新する必要がある
}

void VCFImputable::inverse_haplotype() {
	for(auto p = records.begin(); p != records.end(); ++p) {
		auto	*record = *p;
		record->inverse_haplotype();
	}
}

size_t VCFImputable::max_index() const {
	size_t	index = records.front()->get_index();
	for(auto p = records.begin() + 1; p != records.end(); ++p) {
		if((*p)->get_index() > index)
			index = (*p)->get_index();
	}
	return index;
}

void VCFImputable::which_comes_from_hetero_parent_chromosomes() {
	for(auto p = records.begin(); p != records.end(); ++p)
		(*p)->which_comes_from_hetero_parent_chromosomes();
}

void VCFImputable::walk(size_t v0, vector<vector<int>>& haplo,
									const Graph::InvGraph& tree_graph) {
	auto	q = tree_graph.find(v0);
	for(auto p = q->second.begin(); p != q->second.end(); ++p) {
		const size_t	v = get<0>(*p);
		const bool		b = get<2>(*p);
		if(haplo[0][v] != -1)
			continue;
		if(b) {		// inverse
			haplo[0][v] = haplo[0][v0] == 0 ? 1 : 0;
			haplo[1][v] = haplo[1][v0] == 0 ? 1 : 0;
		}
		else {
			haplo[0][v] = haplo[0][v0];
			haplo[1][v] = haplo[1][v0];
		}
		walk(v, haplo, tree_graph);
	}
}

Graph::InvGraph VCFImputable::make_graph(VCFHeteroHomo *vcf, int max_dist) {
	const size_t	L = vcf->size();
	vector<vector<int>>	gtss;
	for(size_t i = 0U; i < L; ++i) {
		auto	*record = vcf->get_record(i);
		bool mat_hetero = VCFImputable::determine_which_parents_is_hetero(vcf);
		gtss.push_back(record->genotypes_from_hetero_parent(mat_hetero));
	}
	
	vector<double>	cMs(L);
	for(size_t k = 0U; k < L; ++k)
		cMs[k] = vcf->cM(k);
	
	Graph::InvGraph	graph;
	// キーをあらかじめ登録しておかないと、孤立した点が表されない
	for(size_t k = 0U; k < L; ++k)
		graph[k];
	for(size_t k = 0U; k < L; ++k) {
		for(size_t l = k + 1; l < L; ++l) {
			const double	cM = cMs[l] - cMs[k];
			if(cM > 10.0)
				break;
			const auto	p = vcf->distance(gtss[k], gtss[l], max_dist);
			const int	d = p.first;
			const bool	b = p.second;
			if(d <= max_dist) {
				graph[k].push_back(tuple<size_t,int,bool>(l, d, b));
				graph[l].push_back(tuple<size_t,int,bool>(k, d, b));
			}
		}
	}
	return graph;
}

VCFImputable *VCFImputable::make_subvcf(VCFHeteroHomo *vcf,
										const Graph::InvGraph& graph,
										bool is_mat) {
	const auto	haplo_ = make_parent_haplotypes(vcf->size(), graph);
	auto	indices_ = Graph::keys(graph);
	std::sort(indices_.begin(), indices_.end());
	size_t	i = indices_.front();
	size_t	j = indices_.back();
	const double	total_cM = vcf->cM(j) - vcf->cM(i);
	vector<VCFHeteroHomoRecord *>	new_records;
	vector<vector<int>>	haplo(2);
	vector<size_t>	indices;
	vector<bool>	bs;
	for(auto p = indices_.begin(); p != indices_.end(); ++p) {
		auto	*record = vcf->get_record(*p);
		if(record->is_valid_segregation(is_mat, total_cM)) {
			// VCFImputable::createでVCFImputableRecordを作るから
			// ここではcopyしない
			new_records.push_back(record);
			indices.push_back(*p);
			haplo[0].push_back(haplo_[0][p-indices_.begin()]);
			haplo[1].push_back(haplo_[1][p-indices_.begin()]);
		}
	}
	
	if(new_records.empty())
		return NULL;
	
	bool	mat_hetero = VCFImputable::determine_which_parents_is_hetero(vcf);
	// MapはVCFHeteroHomoから引き継がれる
	auto *new_vcf = create(vcf->get_header(), new_records, mat_hetero,
											haplo, indices, vcf->get_map());
	new_vcf->which_comes_from_hetero_parent_chromosomes();
	return new_vcf;
}

// markerを省いたのでindexを付け替える
void VCFImputable::renumber_indices(vector<VCFImputable *>& vcfs) {
	vector<size_t>	indices;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*vcf = *p;
		for(auto q = vcf->records.begin(); q != vcf->records.end(); ++q)
			indices.push_back((*q)->get_index());
	}
	
	map<size_t,size_t>	dic_indices;
	for(size_t i = 0U; i < indices.size(); ++i)
		dic_indices[indices[i]] = i;
	
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		const auto	*vcf = *p;
		for(auto q = vcf->records.begin(); q != vcf->records.end(); ++q) {
			auto	*record = *q;
			record->set_index(dic_indices[record->get_index()]);
		}
	}
}

VCFCollection *VCFImputable::determine_haplotype(VCFHeteroHomo *vcf,
									const OptionImpute *option, bool is_mat) {
	Graph::InvGraph	graph = make_graph(vcf, option->max_dist);
	const auto	subgraphs = Graph::divide_graph_into_connected(graph);
	
	// 小さなgraphのmarkerは省いて、ヘテロ親から0 or 1のどちらが来たかにする
	vector<VCFImputable *>	vcfs;
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		const auto	subgraph = *p;
		if((int)subgraph.size() < option->min_graph)
			continue;
		VCFImputable	*new_vcf = make_subvcf(vcf, subgraph, is_mat);
		if(new_vcf != NULL)
			vcfs.push_back(new_vcf);
	}
	
	renumber_indices(vcfs);
	return new VCFCollection(vcfs);
}

VCFImputable *VCFImputable::create(const vector<STRVEC>& header,
								const vector<VCFHeteroHomoRecord *>& records_,
								bool hetero, const vector<vector<int>>& haplo,
								const vector<size_t>& indices, const Map& m) {
	vector<VCFImputableRecord *>	records;
	vector<int>	empty;
	const STRVEC&	samples = records_.front()->get_samples();
	for(size_t i = 0U; i < records_.size(); ++i) {
		auto	*r = records_[i];
		records.push_back(new VCFImputableRecord(r->get_v(), samples,
												haplo[0][i], haplo[1][i],
												hetero, indices[i], empty));
	}
	return new VCFImputable(header, samples, records, m);
}

bool VCFImputable::determine_which_parents_is_hetero(VCFFamily *vcf) {
	if(vcf->is_all_hetero(true)) {	// is mat hetero?
		if(vcf->is_all_homo(false)) {
			return true;
		}
		else {
			cerr << "error : both parents are hetero." << endl;
			exit(1);
		}
	}
	else if(vcf->is_all_hetero(false)) {
		if(vcf->is_all_homo(true)) {
			return false;
		}
		else {
			cerr << "error : both parents are hetero." << endl;
			exit(1);
		}
	}
	else {
		cerr << "error : no parent is hetero." << endl;
		for(size_t i = 0U; i < vcf->size(); ++i) {
			const auto	*record = vcf->get_record(i);
			const auto&	v = record->get_v();
			cerr << i << " " << v[9] << " " << v[10] << " " << endl;
		}
		exit(1);
	}
}

// Lとgraph.size()が違うときおかしい？
vector<vector<int>> VCFImputable::make_parent_haplotypes(size_t L,
												const Graph::InvGraph& graph) {
	const Graph::InvGraph	tree = Graph::minimum_spanning_tree(graph);
	auto	vs = Graph::keys(tree);
	std::sort(vs.begin(), vs.end());
	vector<vector<int>>	haplo(2, vector<int>(L, -1));
	haplo[0][vs[0]] = 0;
	haplo[1][vs[0]] = 1;
	walk(vs[0], haplo, tree);
	
	vector<vector<int>>	haplo2(2);
	for(int i = 0; i < 2; ++i) {
		for(auto p = haplo[i].begin(); p != haplo[i].end(); ++p) {
			if(*p != -1)
				haplo2[i].push_back(*p);
		}
	}
	return haplo2;
}

void VCFImputable::impute_by_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	auto	vcf = c->vcf;
	const size_t	num = vcf->get_samples().size() - 2;
	for(size_t i = c->first; i < num; i += c->num_thread) {
		const string	seq = vcf->make_seq(i);
		const string	hidden_seq = vcf->impute_each_seq(seq, c->Ts, c->MIN);
		c->hs[i] = hidden_seq;
	}
}
