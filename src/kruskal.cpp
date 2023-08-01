#include <tuple>
#include "../include/kruskal.h"
#include "../include/UnionFind.h"

using namespace std;


//////////////////// Kruskal ////////////////////

vector<Kruskal::Edge> Kruskal::Kruskal_core(const vector<Edge>& edges,
													const GraphBase& graph) {
	const vector<Node>	nodes = graph.collect_nodes();
	UnionFind	tree(nodes);
	
	vector<Edge>	selected_edges;
	size_t	counter = 0U;
	for(auto p = edges.begin(); p != edges.end(); ++p) {
		const Node	v1 = get<0>(*p);
		const Node	v2 = get<1>(*p);
		if(v1 > v2)
			continue;
		const Edge	edge(v1, v2);
		if(tree.root(v1) != tree.root(v2)) {
			tree.join(v1, v2);
			selected_edges.push_back(edge);
			counter += 1;
			if(counter == nodes.size() - 1)
				break;
		}
	}
	
	return selected_edges;
}
