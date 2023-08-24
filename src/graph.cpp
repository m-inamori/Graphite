#include <algorithm>
#include "../include/graph.h"
#include "../include/UnionFind.h"

using namespace std;


//////////////////// GraphBase::EdgeGenerator ////////////////////

GraphBase::EdgeGenerator::EdgeGenerator(const GraphBase& g) : graph(g),
													nodes(g.collect_nodes()) {
	p = nodes.begin();
	if(p == nodes.end())
		return;
	
	neighs = graph.neighbors(*p);
	q = neighs.begin();
	if(!is_effective())
		find_next_edge();
}

bool GraphBase::EdgeGenerator::is_effective() const {
	return q != neighs.end();
}

void GraphBase::EdgeGenerator::proceed() {
	if(q == neighs.end()) {
		++p;
		if(p == nodes.end())
			return;
		neighs = graph.neighbors(*p);
		q = neighs.begin();
	}
	else {
		++q;
	}
}

void GraphBase::EdgeGenerator::find_next_edge() {
	while(!ends()) {
		proceed();
		if(is_effective())
			break;
	}
}

bool GraphBase::EdgeGenerator::ends() const {
	return p == nodes.end();
}

GraphBase::Edge GraphBase::EdgeGenerator::next() {
	Edge	edge(*p, *q);
	find_next_edge();
	return edge;
}

size_t GraphBase::count_groups(const map<GraphBase::Node, size_t>& groups) {
	size_t	max_group_index = 0;
	for(auto p = groups.begin(); p != groups.end(); ++p) {
		if(p->second > max_group_index)
			max_group_index = p->second;
	}
	return max_group_index + 1;
}

vector<vector<size_t>> GraphBase::divide_nodes_into_connected() const {
	const vector<size_t>	nodes = collect_nodes();
	UnionFind	uf(nodes);
	EdgeGenerator	eg(*this);
	while(!eg.ends()) {
		const Edge	edge = eg.next();
		uf.join(edge.first, edge.second);
	}
	
	const map<Node, size_t>	groups = uf.divide_into_trees();
	const size_t	num_groups = count_groups(groups);
	vector<vector<size_t>>	vs(num_groups);
	for(auto p = groups.begin(); p != groups.end(); ++p) {
		const size_t	group_index = p->second;
		vs[group_index].push_back(p->first);
	}
	return vs;
}
