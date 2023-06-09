#include <stack>
#include <algorithm>
#include <cmath>
#include <random>
#include <cassert>
#include "../include/common.h"
#include "../include/VCFHeteroHomo.h"
#include "../include/VCFOriginal.h"
#include "../include/Map.h"
#include "../include/Pedigree.h"
#include "../include/Imputer.h"
#include "../include/inverse_graph.h"
#include "../include/option.h"

using namespace std;


//////////////////// VCFHeteroHomoRecord ////////////////////

FillType VCFHeteroHomoRecord::get_fill_type() const {
	if(this->is_imputable())
		return this->is_mat_hetero() ? FillType::MAT : FillType::PAT;
	else
		return FillType::IMPUTABLE;
}

vector<int> VCFHeteroHomoRecord::genotypes_from_hetero_parent() const {
	const vector<int>	gts = this->get_int_gts();
	const int	homo_gt = this->is_mat_hetero() ? pat_int_gt() : mat_int_gt();
	vector<int>	which_froms;
	for(auto p = gts.begin() + 2; p != gts.end(); ++p) {
		const int	gt = homo_gt == 0 ? *p : *p - 1;
		which_froms.push_back(gt != 0 && gt != 1 ? -1 : gt);
	}
	return which_froms;
}

void VCFHeteroHomoRecord::set_haplo(int h) {
	const int	hetero_col = this->is_mat_hetero() ? 9 : 10;
	const int	homo_col = hetero_col == 9 ? 10 : 9;
	this->v[hetero_col] = h == 0 ? "0|1" : "1|0";
	this->v[homo_col] = this->v[homo_col].c_str()[0] == '0' ? "0|0" : "1|1";
	
	const auto	gts = this->genotypes_from_hetero_parent();
	for(size_t i = 0; i < gts.size(); ++i)
		this->which_comes_from[i] = gts[i] == h ? 0 : 1;
}

void VCFHeteroHomoRecord::set_int_gt_by_which_comes_from(
												const vector<int>& ws) {
	this->which_comes_from = ws;
	const string&	mat_gt = this->v[9];
	const string&	pat_gt = this->v[10];
	if(this->is_mat_hetero()) {
		const string	homo_gt = pat_gt.substr(0, 1);
		for(size_t i = 0; i < this->num_progenies(); ++i) {
			this->set_GT(i+2, mat_gt.substr(ws[i]*2, 1) + '|' + homo_gt);
		}
	}
	else {
		const string	homo_gt = mat_gt.substr(0, 1);
		for(size_t i = 0; i < this->num_progenies(); ++i) {
			this->set_GT(i+2, homo_gt + '|' + pat_gt.substr(ws[i]*2, 1));
		}
	}
}

//////////////////// VCFHeteroHomo ////////////////////

VCFHeteroHomo::VCFHeteroHomo(const vector<STRVEC>& h, const STRVEC& s,
							vector<VCFHeteroHomoRecord *> rs, const Map& m) :
						VCFBase(h, s), VCFSmallBase(),
						VCFFamilyBase(), VCFMeasurable(m), records(rs) { }

VCFHeteroHomo::~VCFHeteroHomo() {
	for(auto p = records.begin(); p != records.end(); ++p)
		delete *p;
}

Graph::InvGraph VCFHeteroHomo::make_graph(double max_dist) const {
	vector<vector<int>>	gtss;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		gtss.push_back((*p)->genotypes_from_hetero_parent());
	}
	
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	Graph::InvGraph	graph;
	// キーをあらかじめ登録しておかないと、孤立した点が表されない
	for(size_t k = 0U; k < L; ++k)
		graph[k];
	for(size_t k = 0; k < L; ++k) {
		for(size_t l = k + 1; l < L; ++l) {

			const double	cM = cMs[l] - cMs[k];
			if(cM > 10.0)	// 10cM以上離れていたら繋がりを見ない
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

// ベータ分布を使うので、distはdouble
pair<double, bool> VCFHeteroHomo::distance(const vector<int>& gts1,
											const vector<int>& gts2) {
	const int	N = (int)gts1.size();
	int	counter_right = 0;
	int	counter_diff = 0;
	int	counter_NA = 0;
	for(size_t i = 0; i < gts1.size(); ++i) {
		const int	int_gt1 = gts1[i];
		const int	int_gt2 = gts2[i];
		if(int_gt1 == -1 || int_gt2 == -1)
			counter_NA += 1;
		else if(int_gt1 == int_gt2)
			counter_right += 1;
		else
			counter_diff += 1;
	}
	
	if(counter_right >= counter_diff)
		return make_pair(dist_with_NA(counter_right, counter_NA, N), false);
	else
		return make_pair(dist_with_NA(counter_diff, counter_NA, N), true);
}

double VCFHeteroHomo::dist_with_NA(int right, int counter_NA, int N) {
	const int	diff = N - right - counter_NA;
	const double	diff_ratio = (diff + 1.0) / (right + diff + 2.0);
	return diff + counter_NA * diff_ratio;
}

// haplotype1の各レコードのGenotypeが0なのか1なのか
vector<int> VCFHeteroHomo::make_parent_haplotypes(
										const Graph::InvGraph& graph) const {
	const Graph::InvGraph	tree = Graph::minimum_spanning_tree(graph);
	auto	vs = Graph::keys(tree);
	std::sort(vs.begin(), vs.end());
	const size_t	L = this->size();
	vector<int>	haplo(L, -1);
	haplo[vs[0]] = 0;
	walk(vs[0], haplo, tree);
	
	vector<int>	haplo2;
	for(auto p = haplo.begin(); p != haplo.end(); ++p) {
		if(*p != -1)
				haplo2.push_back(*p);
	}
	return haplo2;
}

void VCFHeteroHomo::walk(size_t v0, vector<int>& haplo,
									const Graph::InvGraph& tree_graph) {
	// 再帰をやめてstackを明示的に使う
	stack<size_t>	stk;
	stk.push(v0);
	while(!stk.empty()) {
		const size_t	v = stk.top();
		stk.pop();
		const auto	q = tree_graph.find(v);
		const auto&	vs = q->second;
		for(auto p = vs.begin(); p != vs.end(); ++p) {
			const size_t	v1 = get<0>(*p);
			const bool		inv = get<2>(*p);
			if(haplo[v1] != -1)		// visited
				continue;
			if(inv)
				haplo[v1] = haplo[v] == 0 ? 1 : 0;
			else
				haplo[v1] = haplo[v];
			stk.push(v1);
		}
	}
}

VCFHeteroHomo *VCFHeteroHomo::make_subvcf(const Graph::InvGraph& graph) const {
	const vector<int>	haplo = make_parent_haplotypes(graph);
	vector<int>	indices;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		indices.push_back(p->first);
	std::sort(indices.begin(), indices.end());
	vector<VCFHeteroHomoRecord *>	records;
	for(size_t i = 0; i < haplo.size(); ++i) {
		VCFHeteroHomoRecord	*record = this->records[indices[i]];
		record->set_haplo(haplo[i]);
		records.push_back(record);
	}
	return new VCFHeteroHomo(this->header, this->samples,
										records, this->get_map());
}

pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>
		VCFHeteroHomo::determine_haplotype(const OptionImpute *option) const {
	const double	max_dist = std::min((double)option->max_dist,
									this->num_progenies() * 0.1);
	const Graph::InvGraph	graph = this->make_graph(max_dist);
	const vector<Graph::InvGraph>	subgraphs =
									Graph::divide_graph_into_connected(graph);
	
	// 小さなgraphのmarkerは省く
	vector<Graph::InvGraph>	gs;
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		if((int)p->size() >= option->min_graph) {
			gs.push_back(*p);
		}
	}
	
	// 小さなgraphしかなければ、その中でも一番大きなgraphにする
	if(gs.empty()) {
		auto	max_g = std::max_element(subgraphs.begin(), subgraphs.end(),
										[](const auto& a, const auto& b) {
											return a.size() < b.size(); });
		gs.push_back(*max_g);
		for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
			if(p != max_g) {
				for(auto q = p->begin(); q != p->end(); ++q) {
					unused_records.push_back(this->records[q->first]);
				}
			}
		}
	}
	else {
		for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
			if((int)p->size() < option->min_graph) {
				for(auto q = p->begin(); q != p->end(); ++q) {
					unused_records.push_back(this->records[q->first]);
				}
			}
		}
	}
	
	vector<VCFHeteroHomo *>	subvcfs;
	for(auto p = gs.begin(); p != gs.end(); ++p) {
		subvcfs.push_back(this->make_subvcf(*p));
	}
	return make_pair(subvcfs, unused_records);
}

string VCFHeteroHomo::make_seq(size_t i) const {
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

bool VCFHeteroHomo::is_all_same_without_N(const string& seq) {
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

string VCFHeteroHomo::create_same_color_string(const string& seq) {
	char	c = '0';	// dummy
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		if(*p != 'N') {
			c = *p;
			break;
		}
	}
	return std::string(seq.size(), c);
}

string VCFHeteroHomo::impute_each_sample_seq(int i,
								const vector<double>& cMs, double min_c) {
	const string	seq = this->make_seq(i);
	if(is_all_same_without_N(seq))
		return create_same_color_string(seq);
	
	const vector<char>	hidden_states = { '0', '1' };
	const vector<char>	states = Imputer::create_states(seq);
	const string	hidden_seq = Imputer::impute(seq,
												hidden_states, states, cMs);
	const string	painted_seq = Imputer::paint(hidden_seq, cMs, min_c);
	return painted_seq;
}

void VCFHeteroHomo::impute_each(const OptionImpute *option) {
	vector<double>	cMs;
	const size_t	L = this->size();
	for(size_t i = 0; i < L; ++i) {
		cMs.push_back(this->record_cM(i));
	}
	
	vector<string>	imputed_seqs;
	for(size_t i = 0; i < this->num_progenies(); ++i) {
		imputed_seqs.push_back(
					impute_each_sample_seq(i, cMs, option->min_crossover));
	}
	
	for(size_t k = 0; k < this->size(); ++k) {
		VCFHeteroHomoRecord	*record = this->records[k];
		vector<int>	ws;
		for(auto p = imputed_seqs.begin(); p != imputed_seqs.end(); ++p)
			ws.push_back(p->c_str()[k] - '0');
		record->set_int_gt_by_which_comes_from(ws);
	}
}

const OptionImpute *VCFHeteroHomo::create_option() const {
	const size_t	num = std::max(2UL, this->size());
	const double	cM_length = this->record_cM(records.size()-1);
	const int	max_dist = std::max(4,
						(int)(cM_length * num_progenies()
										/ num * log10(num) * 2.5 * 0.01));
	return new OptionImpute(max_dist, 20, 5, 1.0);
}

pair<vector<VCFHeteroHomo *>, vector<VCFHeteroHomoRecord *>>
													VCFHeteroHomo::impute() {
	if(this->records.empty()) {
		// vcfを新たに作らないと、deleteしていけないvcfをdeleteしてしまう
		VCFHeteroHomo	*empty_vcf = new VCFHeteroHomo(header, samples,
								vector<VCFHeteroHomoRecord *>(), get_map());
		return make_pair(vector<VCFHeteroHomo *>(1, empty_vcf),
									vector<VCFHeteroHomoRecord *>());
	}
	
	const OptionImpute	*option = this->create_option();
	auto	p = this->determine_haplotype(option);
	vector<VCFHeteroHomo *>&	vcfs = p.first;
	vector<VCFHeteroHomoRecord *>&	unused_records = p.second;
	for(auto q = vcfs.begin(); q != vcfs.end(); ++q)
		(*q)->impute_each(option);
	delete option;
	return make_pair(vcfs, unused_records);
}

// 共通のヘテロ親はどれだけマッチしているか
pair<int, int> VCFHeteroHomo::match(const VCFHeteroHomo *other) const {
	if(this->records.empty() || other->records.empty())
		return pair<int, int>(0, 0);
	
	const string	mat_GT1 = this->records.front()->get_GT(0);
	const string	mat_GT2 = other->records.front()->get_GT(0);
	const int	hetero_col1 = mat_GT1.c_str()[0] != mat_GT1.c_str()[2] ? 9 : 10;
	const int	hetero_col2 = mat_GT2.c_str()[0] != mat_GT2.c_str()[2] ? 9 : 10;
	int	num_match = 0;
	int	num_unmatch = 0;
	size_t	k = 0;
	size_t	l = 0;
	while(k < this->size() && l < other->size()) {
		const VCFHeteroHomoRecord	*record1 = this->records[k];
		const VCFHeteroHomoRecord	*record2 = other->records[l];
		if(record1->pos() == record2->pos()) {
			if(record1->get_v()[hetero_col1] == record2->get_v()[hetero_col2])
				num_match += 1;
			else
				num_unmatch += 1;
		}
		if(record1->pos() <= record2->pos())
			k += 1;
		if(record1->pos() >= record2->pos())
			l += 1;
	}
	
	return make_pair(num_match, num_unmatch);
}

void VCFHeteroHomo::inverse_hetero_parent_phases() {
	const int	hetero_index = this->is_mat_hetero() ? 0 : 1;
	for(auto p = this->records.begin(); p != this->records.end(); ++p) {
		if((*p)->get_GT(hetero_index) == "0|1")
			(*p)->set_GT(hetero_index, "1|0");
		else
			(*p)->set_GT(hetero_index, "0|1");
	}
}

void VCFHeteroHomo::impute_in_thread(void *config) {
	const auto	*c = (ConfigThread *)config;
	const size_t	n = c->size();
	// 家系数が少ない場合、スレッドの内部でも並列に処理したい
	const int	num_inner_threads = (c->option->num_threads + n - 1) / n;
	for(size_t i = c->first; i < n; i += c->num_threads) {
		auto	result = impute_vcfs(c->vcfs_heho[i],
											c->option, num_inner_threads);
		c->imputed_vcfs[i] = result;
	}
}

vector<ImpResult> VCFHeteroHomo::impute_hetero_homo_all(
						const map<string, vector<VCFHeteroHomo *>>& vcfs,
						const Option *option) {
	vector<ImpResult>	results(vcfs.size());
	
	// vectorにしたほうがマルチスレッドにしやすい
	vector<vector<VCFHeteroHomo *>>	vcfss;
	for(auto p = vcfs.begin(); p != vcfs.end(); ++p) {
		vcfss.push_back(p->second);
	}
	
	const int	T = min((int)vcfs.size(), option->num_threads);
	vector<ConfigThread *>	configs(T);
	for(int i = 0; i < T; ++i)
		configs[i] = new ConfigThread(vcfss, option, i, T, results);
	
#ifndef DEBUG
	vector<pthread_t>	threads_t(T);
	for(int i = 0; i < T; ++i)
		pthread_create(&threads_t[i], NULL,
					(void *(*)(void *))&impute_in_thread, (void *)configs[i]);
	
	for(int i = 0; i < T; ++i)
		pthread_join(threads_t[i], NULL);
#else
	for(int i = 0; i < T; ++i)
		impute_in_thread(configs[i]);
#endif
	
	Common::delete_all(configs);
	
	return results;
}


// ヘテロ親が同じVCFを集めて補完する
// ついでにphaseもなるべく同じになるように変更する
ImpResult VCFHeteroHomo::impute_vcfs(const vector<VCFHeteroHomo *>& vcfs,
										const Option* op, int num_threads) {
	vector<VCFHeteroHomo *>	imputed_vcfs;
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto q = vcfs.begin(); q != vcfs.end(); ++q) {
		auto	*vcf = *q;
		auto	p = vcf->impute();
		auto	subvcfs = p.first;
		auto	unused = p.second;
		imputed_vcfs.insert(imputed_vcfs.end(), subvcfs.begin(), subvcfs.end());
		unused_records.insert(unused_records.end(),
										unused.begin(), unused.end());
		// 元のVCFは消したいが、Recordsは使いまわししているので、
		// 元のVCFのRecordsを消してからVCFのdeleteをする
		vcf->records.clear();
		delete vcf;
	}
	inverse_phases(imputed_vcfs);
	return make_pair(imputed_vcfs, unused_records);
}

void VCFHeteroHomo::inverse_phases(const vector<VCFHeteroHomo *>& vcfs) {
	
	const auto	*graph = make_vcf_graph(vcfs);
	const vector<bool>	invs = graph->optimize_inversions();
	delete graph;
	for(size_t i = 0; i < vcfs.size(); ++i) {
		if(invs[i])
			vcfs[i]->inverse_hetero_parent_phases();
	}
}

const InverseGraph *VCFHeteroHomo::make_vcf_graph(
									const vector<VCFHeteroHomo *>& vcfs) {
	const size_t	N = vcfs.size();
	vector<vector<tuple<size_t, int, int>>>	graph(N);
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = i + 1; j < N; ++j) {
			const auto	p = vcfs[i]->match(vcfs[j]);
			const int	num_match = p.first;
			const int	num_unmatch = p.second;
			if(num_match != 0 || num_unmatch != 0) {
				graph[i].push_back(make_tuple(j, num_match, num_unmatch));
				graph[j].push_back(make_tuple(i, num_match, num_unmatch));
			}
		}
	}
	return InverseGraph::convert(graph);
}

vector<vector<bool>> VCFHeteroHomo::enumerate_bools(size_t L) {
	vector<vector<bool>>	bss;
	if(L <= 10) {
		for(size_t n = 0; n < (1U << L); ++n) {
			vector<bool>	bs(L);
			for(size_t k = 0; k < L; ++k)
				bs[k] = ((n >> k) & 1) == 1;
			bss.push_back(bs);
		}
	}
	else {
		std::mt19937	mt(123456789);
		std::uniform_int_distribution<int>	dist(0, 1);
		
		for(size_t i = 0; i < 1024; ++i) {
			vector<bool>	bs(L);
			for(size_t j = 1; j < L; ++j)
				bs[j] = dist(mt) == 0;
			bss.push_back(bs);
		}
	}
	return bss;
}

int VCFHeteroHomo::match_score(const vector<bool>& invs,
					const vector<vector<tuple<int, int, int>>>& graph) {
	const size_t	N = graph.size();
	int	num_wrong = 0;
	for(int i = 0; i < (int)N; ++i) {
		const auto&	v = graph[i];
		for(auto q = v.begin(); q != v.end(); ++q) {
			const int	j = std::get<0>(*q);
			const int	n1 = std::get<1>(*q);
			const int	n2 = std::get<2>(*q);
			if(i > j)
				continue;
			const bool	inv1 = i == 0 ? false : invs[i-1];
			const bool	inv2 = invs[j-1];
			if(inv1 == inv2)
				num_wrong += n2;
			else
				num_wrong += n1;
		}
	}
	return num_wrong;
}

// FamilyごとにVCFHeteroHomoを作って親ごとに格納する
tuple<VCFHeteroHomo *, VCFHeteroHomo *, vector<VCFHeteroHomoRecord *>>
VCFHeteroHomo::make_VCFHeteroHomo(const vector<VCFHeteroHomoRecord *>& records,
								  const Family *family,
								  const VCFSmall *vcf, const Map& geno_map) {
	const STRVEC&	samples = family->get_samples();
	const vector<STRVEC>	header = vcf->create_header(samples);
	vector<VCFHeteroHomoRecord *>	heho_mat_records;
	vector<VCFHeteroHomoRecord *>	heho_pat_records;
	vector<VCFHeteroHomoRecord *>	unused_records;
	for(auto p = records.begin(); p != records.end(); ++p) {
		VCFHeteroHomoRecord	*record = *p;
		if(!record->is_imputable())
			unused_records.push_back(record);
		else if(record->is_mat_hetero())
			heho_mat_records.push_back(record);
		else
			heho_pat_records.push_back(record);
	}
	auto	*vcf_mat = new VCFHeteroHomo(header, samples,
											heho_mat_records, geno_map);
	auto	*vcf_pat = new VCFHeteroHomo(header, samples,
											heho_pat_records, geno_map);
	return make_tuple(vcf_mat, vcf_pat, unused_records);
}

