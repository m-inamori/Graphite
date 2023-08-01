#ifndef __KRUSKAL
#define __KRUSKAL

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "graph.h"


namespace Kruskal {
	
	using Node = std::size_t;
	using Edge = std::pair<Node, Node>;
	
	std::vector<Edge> Kruskal_core(const std::vector<Edge>& edges,
												const GraphBase& graph);
	
	template<typename T>
	std::vector<Edge> Kruskal(const WeightedGraphBase<T>& graph) {
		auto	weighted_edges = graph.collect_weighted_edges();
		std::stable_sort(weighted_edges.begin(), weighted_edges.end(),
								typename WeightedGraphBase<T>::LessWeight());
		std::vector<Edge>	edges;
		for(auto p = weighted_edges.begin(); p != weighted_edges.end(); ++p)
			edges.push_back(Edge(std::get<0>(*p), std::get<1>(*p)));
		
		return Kruskal_core(edges, graph);
	}
}
#endif
