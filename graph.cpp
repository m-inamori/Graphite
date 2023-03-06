#include <iostream>
#include <cassert>
#include "graph.h"
#include "kruskal.h"

using namespace std;


vector<size_t> Graph::keys(const InvGraph& graph) {
	vector<size_t>	vs;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		vs.push_back(p->first);
	return vs;
}

int Graph::get_edge_value(size_t v1, size_t v2, const InvGraph& graph) {
	auto	p = graph.find(v1);
	if(p == graph.end())
		return -1;
	
	const auto&	vs = p->second;
	for(auto q = vs.begin(); q != vs.end(); ++q) {
		if(get<0>(*q) == v2)
			return get<1>(*q);
	}
	return -1;
}

Graph::WeightedGraph Graph::trim_inverse(const InvGraph& graph) {
	WeightedGraph	weighted_graph;
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		weighted_graph[p->first];
	}
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const auto&	v = p->second;
		for(auto q = v.begin(); q != v.end(); ++q)
			weighted_graph[p->first].push_back(
								pair<size_t,int>(get<0>(*q), get<1>(*q)));
	}
	return weighted_graph;
}

Graph::InvGraph Graph::filter_graph(const InvGraph& graph,
									const WeightedGraph& tree) {
	set<pair<size_t,size_t>>	edges;
	for(auto p = tree.begin(); p != tree.end(); ++p) {
		const size_t	v1 = p->first;
		const auto&	vs = p->second;
		for(auto q = vs.begin(); q != vs.end(); ++q)
			edges.insert(pair<size_t,size_t>(v1, q->first));
	}
	
	InvGraph	tree_graph;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		tree_graph[p->first];
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const size_t	v1 = p->first;
		const auto	vs = p->second;
		for(auto q = vs.begin(); q != vs.end(); ++q) {
			const size_t	v2 = get<0>(*q);
			if(edges.find(pair<size_t,size_t>(v1, v2)) != edges.end())
				tree_graph[v1].push_back(*q);
		}
	}
	return tree_graph;
}

Graph::InvGraph Graph::minimum_spanning_tree(const InvGraph& graph) {
	const WeightedGraph	weighted_graph = trim_inverse(graph);
	const WeightedGraph	tree = Kruskal::Kruskal(weighted_graph);
	return filter_graph(graph, tree);
}

void Graph::walk(size_t v0, const InvGraph& graph,
						vector<size_t>& vs, set<size_t>& visited) {
	vs.push_back(v0);
	visited.insert(v0);
	auto	iter = graph.find(v0);
	assert(iter != graph.end());
	const vector<std::tuple<std::size_t,double,bool>>&	vec = iter->second;
	for(auto p = vec.begin(); p != vec.end(); ++p) {
		const size_t&	v = get<0>(*p);
		if(visited.find(v) != visited.end())
			continue;
		walk(v, graph, vs, visited);
	}
}

vector<Graph::InvGraph> Graph::divide_graph_into_connected(
											const InvGraph& graph) {
	vector<Graph::InvGraph>	graphs;
	set<size_t>	visited;
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const size_t&	v = p->first;
		if(visited.find(v) != visited.end())
			continue;
		vector<size_t>	vs;
		walk(v, graph, vs, visited);
		Graph::InvGraph	g;
		for(auto q = vs.begin(); q != vs.end(); ++q) {
			const size_t&	v1 = *q;
			auto	iter = graph.find(v1);
			assert(iter != graph.end());
			g[v1] = iter->second;
		}
		graphs.push_back(g);
	}
	return graphs;
}

// 頭より末尾にpushする方が速いので、逆向きのpathを返す
vector<size_t> Graph::walk(size_t v1, size_t v2, const InvGraph& graph,
													set<size_t>& visited) {
	if(v1 == v2)
		return vector<size_t>(1, v2);
	
	visited.insert(v1);
	auto	p = graph.find(v1);
	if(p == graph.end()) {
		cout << "error : " << v1 << " " << v2 << endl;
		exit(1);
	}
	const auto&	vs = p->second;
	for(auto q = vs.begin(); q != vs.end(); ++q) {
		const size_t	v = get<0>(*q);
		if(visited.find(v) != visited.end())
			continue;
		auto	path = walk(v, v2, graph, visited);
		if(!path.empty()) {		// v2 found
			path.push_back(v1);
			return path;
		}
	}
	visited.erase(v1);
	return vector<size_t>();	// v2 not found
}

vector<size_t> Graph::find_path(size_t v1, size_t v2, const InvGraph& graph) {
	set<size_t>	visited;
	const auto	rev_path = walk(v1, v2, graph, visited);
	return vector<size_t>(rev_path.rbegin(), rev_path.rend());
}
