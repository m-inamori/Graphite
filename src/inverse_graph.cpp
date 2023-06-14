#include <set>
#include <stack>
#include <algorithm>
#include <cassert>
#include "../include/inverse_graph.h"
#include "../include/random.h"

using namespace std;


vector<InverseGraph::Node> InverseGraph::sort_nodes() const {
	vector<Node>	vs;
	for(auto p = g.begin(); p != g.end(); ++p)
		vs.push_back(p->first);
	std::sort(vs.begin(), vs.end());
	return vs;
}

map<InverseGraph::Node, size_t> InverseGraph::calc_dic() const {
	map<Node, size_t>	dic;
	for(size_t i = 0; i < vs.size(); ++i)
		dic[vs[i]] = i;
	return dic;
}

void InverseGraph::walk(size_t v0, vector<size_t>& vs,
									set<size_t>& visited) const {
	vs.push_back(v0);
	visited.insert(v0);
	auto	iter = g.find(v0);
	assert(iter != g.end());
	const vector<std::tuple<Node, int, int>>&	vec = iter->second;
	for(auto p = vec.begin(); p != vec.end(); ++p) {
		const size_t&	v = get<0>(*p);
		if(visited.find(v) != visited.end())
			continue;
		this->walk(v, vs, visited);
	}
}

vector<bool> InverseGraph::optimize_inversions() const {
	const auto	subgraphs = this->divide_graph_into_connected();
	vector<bool>	bs(this->size(), false);
	for(auto p = subgraphs.begin(); p != subgraphs.end(); ++p) {
		const InverseGraph	*subg = *p;
		const auto	dic_bs = subg->optimize_inversions_connected();
		for(auto q = dic_bs.begin(); q != dic_bs.end(); ++q)
			bs[q->first] = q->second;
		delete subg;
	}
	return bs;
}

// 本当はgraph.cppと共通にしたい
// weightをメソッドにすればよい
vector<const InverseGraph *> InverseGraph::divide_graph_into_connected() const {
	vector<const InverseGraph *>	graphs;
	set<size_t>	visited;
	for(auto p = g.begin(); p != g.end(); ++p) {
		const size_t&	v = p->first;
		if(visited.find(v) != visited.end())
			continue;
		Graph	g1;
		vector<size_t>	vs;
		this->walk(v, vs, visited);
		for(auto q = vs.begin(); q != vs.end(); ++q) {
			const Node&	v1 = *q;
			auto	iter = g.find(v1);
			assert(iter != g.end());
			g1[v1] = iter->second;
		}
		graphs.push_back(new InverseGraph(g1));
	}
	return graphs;
}

// thisは連結成分である前提
map<InverseGraph::Node, bool>
					InverseGraph::optimize_inversions_connected() const {
	if(this->size() <= 20)
		return this->search_all();
	
	// trueとfalseの偏りが大きいものからedgeを確定させていく
	// グラフが全て繋がるまでに矛盾が無ければOK
	const auto	dic_bs = this->connect_biased_edges();
	if(!dic_bs.empty())
		return dic_bs;
	
	// それでもダメならランダム
	// 本当はSAを使いたい
	return this->search_randomly();
}

// しらみつぶし
map<InverseGraph::Node, bool> InverseGraph::search_all() const {
	const size_t	N = this->size();
	int	min_score = 100000000;
	vector<bool>	min_invs;
	for(int i = 0; i < (1 << (N-1)); ++i) {
		// 本当は全てintのままで扱えれば速い
		vector<bool>	invs(N-1);
		for(size_t j = 0; j < N - 1; ++j)
			invs[j] = ((i >> j) & 1) == 1;
		
		const int	score = this->calc_match_score(invs);
		if(score < min_score) {
			min_score = score;
			min_invs = invs;
		}
	}
	
	map<Node, bool>	dic;
	dic[vs[0]] = false;
	for(size_t i = 1; i < vs.size(); ++i)
		dic[vs[i]] = min_invs[i-1];
	return dic;
}

int InverseGraph::choice_seed() const {
	int	s = 0;
	int	D = 1000000000;
	for(auto p = g.begin(); p != g.end(); ++p) {
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			s = (s + (int)(get<0>(*q))) + D;
	}
	return s;
}

std::map<InverseGraph::Node, bool> InverseGraph::search_randomly() const {
	vector<bool>	min_bs;
	int				min_score = 1000000000;
	Random	random(this->choice_seed());
	for(int i = 0; i < 524288; ++i) {	// 2^19
		vector<bool>	bs(this->size() - 1);
		for(size_t k = 0; k < this->size() - 1; ++k)
			bs.push_back(random.randbool());
		const int	cur_score = this->calc_match_score(bs);
		if(min_bs.empty() || cur_score < min_score) {
			min_bs = bs;
			min_score = cur_score;
		}
	}
	
	map<Node, bool>	dic;
	dic[vs[0]] = false;
	for(size_t i = 1; i < vs.size(); ++i)
		dic[vs[i]] = min_bs[i-1];
	return dic;
}

vector<InverseGraph::Edge> InverseGraph::extract_edges() const {
	vector<Edge>	edges;
	for(auto p = g.begin(); p != g.end(); ++p) {
		const Node&	v1 = p->first;
		for(auto q = p->second.begin(); q != p->second.end(); ++q) {
			const Node&	v2 = get<0>(*q);
			const int&	n1 = get<1>(*q);
			const int&	n2 = get<2>(*q);
			if(v1 < v2)
				edges.push_back(Edge(v1, v2, n1, n2));
		}
	}
	return edges;
}

int InverseGraph::calc_match_score(const vector<bool>& bs) const {
	int	num_wrong = 0;
	const auto	edges = this->extract_edges();
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		auto	p1 = dic_vs.find(get<0>(*p));
		auto	p2 = dic_vs.find(get<1>(*p));
		const size_t	i = p1->second;
		const size_t	j = p2->second;
		const bool	b1 = i == 0 ? false : bs[i-1];
		const bool	b2 = bs[j-1];
		if(b1 == b2)
			num_wrong += get<3>(*p);
		else
			num_wrong += get<2>(*p);
	}
	return num_wrong;
}

vector<pair<double, InverseGraph::Edge>> InverseGraph::sort_edges() const {
	const vector<Edge>	edges = this->extract_edges();
	vector<pair<double, Edge>>	sorted_edges;
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const int	n1 = get<2>(*p);
		const int	n2 = get<3>(*p);
		const double	r = (n1 + 1) / (double)(n1 + n2 + 2);
		const double	ratio = std::min(r, 1.0 - r);
		sorted_edges.push_back(make_pair(ratio, *p));
	}
	std::sort(sorted_edges.begin(), sorted_edges.end());
	return sorted_edges;
}

void InverseGraph::join(BoolGraph& graph1, const BoolGraph& graph2) {
	for(auto p = graph2.begin(); p != graph2.end(); ++p)
		graph1.insert(*p);
}

map<InverseGraph::Node, bool> InverseGraph::connect_biased_edges() const {
	const vector<pair<double, Edge>>	sorted_edges = this->sort_edges();
	// [ { node: [(index, inv?)] } ]
	vector<map<size_t, vector<pair<size_t, bool>>>>	subgraphs;
	for(auto p = g.begin(); p != g.end(); ++p) {
		map<size_t, vector<pair<size_t, bool>>>	m;
		m[p->first];
		subgraphs.push_back(m);
	}
	
	// 偏りが大きいエッジから追加していく
	for(auto p = sorted_edges.begin(); p != sorted_edges.end(); ++p) {
		const size_t	i = get<0>(p->second);
		const size_t	j = get<1>(p->second);
		const int		n1 = get<2>(p->second);
		const int		n2 = get<3>(p->second);
		int	k1;
		for(k1 = 0; ; ++k1) {
			auto	g1 = subgraphs[k1];
			if(g1.find(i) != g1.end())
				break;
		}
		int	k2;
		for(k2 = 0; ; ++k2) {
			auto	g2 = subgraphs[k2];
			if(g2.find(j) != g2.end())
				break;
		}
		if(k1 == k2) {
			auto	subg1 = subgraphs[k1];
			subg1[i].push_back(make_pair(j, n1 < n2));
			subg1[j].push_back(make_pair(i, n1 < n2));
			if(!InverseGraph::is_consistent(subgraphs[k1]))
				return map<size_t, bool>();
		}
		else {
			auto&	subg1 = subgraphs[k1];
			auto&	subg2 = subgraphs[k2];
			join(subg1, subg2);
			subg1[i].push_back(make_pair(j, n1 < n2));
			subg1[j].push_back(make_pair(i, n1 < n2));
			subgraphs.erase(subgraphs.begin() + k2);
			if(subgraphs.size() == 1)
				return InverseGraph::invs(subgraphs[0]);
		}
	}
	return map<size_t, bool>();		// ここには来ないはず
}

const InverseGraph *InverseGraph::convert(
						const vector<vector<tuple<Node, int, int>>>& graph) {
	Graph	g;
	for(size_t i = 0; i < graph.size(); ++i) {
		g[i] = graph[i];
	}
	return new InverseGraph(g);
}

bool InverseGraph::is_consistent(
						const map<Node, vector<pair<Node, bool>>>& graph) {
	map<Node, bool>	visited;
	Node	v0 = graph.begin()->first;
	stack<pair<Node, bool>>	stk;
	stk.push(make_pair(v0, false));
	while(!stk.empty()) {
		const Node	v = stk.top().first;
		const bool	b = stk.top().second;
		stk.pop();
		auto	iter = graph.find(v);
		for(auto p = iter->second.begin(); p != iter->second.end(); ++p) {
			const Node	v1 = p->first;
			const bool	b1_edge = p->second;
			const bool	b1 = b ^ b1_edge;
			auto	q = visited.find(v1);
			if(q != visited.end()) {
				const bool	b1_prev = q->second;
				if(b1 != b1_prev)
					return false;
			}
			else {
				visited.insert(make_pair(v1, b1));
				stk.push(make_pair(v1, b1));
			}
		}
	}
	return true;
}

map<InverseGraph::Node, bool> InverseGraph::invs(
							const map<Node, vector<pair<Node, bool>>>& graph) {
	map<Node, bool>	visited;
	const Node	v0 = graph.begin()->first;
	stack<pair<Node, bool>>	stk;
	stk.push(make_pair(v0, false));		// v0は反転しないとする
	while(!stk.empty()) {
		const Node	v = stk.top().first;
		const bool	b = stk.top().second;
		stk.pop();
		auto	iter = graph.find(v);
		for(auto p = iter->second.begin(); p != iter->second.end(); ++p) {
			const Node	v1 = p->first;
			const bool	b1_edge = p->second;
			const bool	b1 = b ^ b1_edge;
			if(visited.find(v1) == visited.end()) {
				visited.insert(make_pair(v1, b1));
				stk.push(make_pair(v1, b1));
			}
		}
	}
	return visited;
}
