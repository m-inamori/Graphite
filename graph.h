#ifndef __GRAPH
#define __GRAPH

#include <vector>
#include <map>
#include <set>

namespace Graph {

typedef std::map<std::size_t,std::vector<std::pair<std::size_t,int>>>
															WeightedGraph;
// graph with weight and inverse
typedef std::map<std::size_t,std::vector<std::tuple<std::size_t,int,bool>>>
																	InvGraph;

std::vector<std::size_t> keys(const InvGraph& graph);

WeightedGraph trim_inverse(const InvGraph& graph);
InvGraph filter_graph(const InvGraph& graph, const WeightedGraph& tree);
InvGraph minimum_spanning_tree(const InvGraph& graph);
void walk(std::size_t v0, const InvGraph& graph,
				std::vector<std::size_t>& vs, std::set<std::size_t>& visited);
std::vector<InvGraph> divide_graph_into_connected(const InvGraph& graph);

}
#endif
