#include <stack>
#include <set>
#include "../include/invgraph.h"
#include "../include/kruskal.h"

using namespace std;

vector<InvGraph::Node> InvGraph::collect_nodes() const {
	vector<Node>	nodes;
	for(auto p = begin(); p != end(); ++p)
		nodes.push_back(p->first);
	return nodes;
}

vector<InvGraph::Node> InvGraph::neighbors(Node v0) const {
	auto	p = find(v0);
	if(p == end())
		return vector<Node>();
	else {
		vector<Node>	neighs;
		for(auto q = p->second.begin(); q != p->second.end(); ++q)
			neighs.push_back(get<0>(*q));
		return neighs;
	}
}

vector<WeightedGraphBase<double>::Edge>
					InvGraph::collect_weighted_edges() const {
	vector<WeightedGraphBase<double>::Edge>	edges;
	for(auto p = begin(); p != end(); ++p) {
		const Node	u = p->first;
		const auto&	neis = p->second;
		for(auto q = neis.begin(); q != neis.end(); ++q) {
			const Node		v = get<0>(*q);
			const double	w = get<1>(*q);
			if(u < v)
				edges.push_back(make_tuple(u, v, w));
		}
	}
	return edges;
}

InvGraph InvGraph::minimum_spanning_tree() const {
	map<pair<Node, Node>, tuple<Node, double, bool>>	m;
	for(auto p = begin(); p != end(); ++p) {
		const Node	u = p->first;
		const auto&	neis = p->second;
		for(auto q = neis.begin(); q != neis.end(); ++q) {
			const Node	v = get<0>(*q);
			m[make_pair(u, v)] = *q;
		}
	}
	
	InvGraph	graph;
	for(auto p = begin(); p != end(); ++p) {
		const Node	u = p->first;
		graph[u];
	}
	
	const auto	edges = Kruskal::Kruskal(*this);
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const auto&	edge = *p;
		auto	q = m.find(edge);
		const Node		u = edge.first;
		const Node		v = edge.second;
		const double	w = get<1>(q->second);
		const bool		b = get<2>(q->second);
		graph[u].push_back(make_tuple(v, w, b));
		graph[v].push_back(make_tuple(u, w, b));
	}
	return graph;
}

vector<InvGraph> InvGraph::divide_into_connected() const {
	const vector<vector<Node>>	nss = divide_nodes_into_connected();
	vector<InvGraph>	graphs(nss.size());
	for(size_t i = 0; i < nss.size(); ++i) {
		for(auto p = nss[i].begin(); p != nss[i].end(); ++p) {
			const Node	v = *p;
			auto	q = find(v);
			graphs[i][v] = q->second;
		}
	}
	return graphs;
}

vector<InvGraph::Edge> InvGraph::walk(size_t v0) const {
	vector<Edge>	edges;
	set<Node>		visited;
	stack<size_t>	stk;
	stk.push(v0);
	visited.insert(v0);
	while(!stk.empty()) {
		const size_t	v = stk.top();
		stk.pop();
		const auto	q = find(v);
		const auto&	vs = q->second;
		for(auto p = vs.begin(); p != vs.end(); ++p) {
			const Node	v1 = get<0>(*p);
			const bool	inv = get<2>(*p);
			const Edge	edge = make_tuple(v, v1, inv);
			edges.push_back(edge);
			if(visited.find(v1) != visited.end())
				continue;
			stk.push(v1);
			visited.insert(v1);
		}
	}
	return edges;
}
