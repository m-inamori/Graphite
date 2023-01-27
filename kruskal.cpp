#include <tuple>
#include <algorithm>
#include "kruskal.h"

using namespace std;


//////////////////// UnionFind ////////////////////

Kruskal::UnionFind::UnionFind(const vector<size_t>& nodes) {
	for(auto p = nodes.begin(); p != nodes.end(); ++p) {
		parents.insert(pair<size_t,size_t>(*p, *p));
		heights.insert(pair<size_t,int>(*p, 1));
	}
}

void Kruskal::UnionFind::join(const size_t& v1, const size_t& v2) {
	const size_t&	r1 = this->root(v1);
	const size_t&	r2 = this->root(v2);
	const int	h1 = this->heights[r1];
	const int	h2 = this->heights[r2];
	if(h1 <= h2) {
		this->parents[r1] = r2;
		this->heights[r2] = std::max(this->heights[r2], this->heights[r1]+1);
	}
	else {
		this->parents[r2] = r1;
		this->heights[r1] = std::max(this->heights[r1], this->heights[r2]+1);
	}
}

const size_t& Kruskal::UnionFind::root(const size_t& v0) const {
	size_t	v = v0;
	while(true) {
		auto	p = this->parents.find(v);
		if(p->second == v)
			return p->second;
		v = p->second;
	}
}

//////////////////// Kruskal ////////////////////

Kruskal::GRAPH Kruskal::Kruskal(const GRAPH& graph) {
	vector<size_t>	nodes;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		nodes.push_back(p->first);
	
	UnionFind	tree(nodes);
	
	vector<tuple<size_t,size_t,int>>	edges;
	for(auto p = graph.begin(); p != graph.end(); ++p) {
		const size_t&	v1 = p->first;
		const vector<pair<size_t,double>>&	vec = p->second;
		for(auto q = vec.begin(); q != vec.end(); ++q) {
			const size_t&	v2 = q->first;
			const int		w = q->second;
			if(v1 < v2)
				edges.push_back(tuple<size_t,size_t,int>(v1, v2, w));
		}
	}
	std::stable_sort(edges.begin(), edges.end(),
					[](const tuple<size_t,size_t,int>& a,
					   const tuple<size_t,size_t,int>& b) {
								return get<2>(a) < get<2>(b);
						});
	
	GRAPH	new_graph;
	for(auto p = graph.begin(); p != graph.end(); ++p)
		new_graph[p->first];
	size_t	counter = 0U;
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const size_t&	v1 = get<0>(*p);
		const size_t&	v2 = get<1>(*p);
		double			w  = get<2>(*p);
		if(tree.root(v1) != tree.root(v2)) {
			tree.join(v1, v2);
			new_graph[v1].push_back(pair<size_t,int>(v2, w));
			new_graph[v2].push_back(pair<size_t,int>(v1, w));
			counter += 1;
			if(counter == nodes.size() - 1)
				break;
		}
	}
	
	return new_graph;
}
