#ifndef __INVERSEGRAPH
#define __INVERSEGRAPH

#include <vector>
#include <map>
#include <set>
#include <tuple>


//////////////////// InverseGraph ////////////////////

class InverseGraph {
public:
	typedef std::size_t	Node;
	typedef std::map<Node, std::vector<std::tuple<Node, int, int>>>	Graph;
	typedef std::tuple<Node, Node, int, int>	Edge;
	typedef std::map<Node, std::vector<std::pair<Node, bool>>>	BoolGraph;
	
private:
	const Graph	g;
	const std::vector<Node>	vs;		// Nodeを昇順に並べた配列
	const std::map<Node, std::size_t>	dic_vs;		// Node -> index
	
public:
	InverseGraph(const Graph& g_) : g(g_), vs(sort_nodes()),
											dic_vs(calc_dic()) { }
	
	std::size_t	size() const { return this->g.size(); }
	
	std::vector<bool> optimize_inversions() const;
	std::vector<const InverseGraph *> divide_graph_into_connected() const;
	std::map<Node, bool> optimize_inversions_connected() const;
	
private:
	std::vector<Node> sort_nodes() const;
	std::map<Node, std::size_t> calc_dic() const;
	void walk(Node v0, std::vector<Node>& vs,
										std::set<Node>& visited) const;
	std::map<Node, bool> search_all() const;
	int choice_seed() const;
	std::map<Node, bool> search_randomly() const;
	std::vector<Edge> extract_edges() const;
	int calc_match_score(const std::vector<bool>& bs) const;
	std::vector<std::pair<double, Edge>> sort_edges() const;
	std::map<Node, bool> connect_biased_edges() const;
	
public:
	static const InverseGraph *convert(
	const std::vector<std::vector<std::tuple<Node, int, int>>>& graph);
	static bool is_consistent(
			const std::map<Node, std::vector<std::pair<Node, bool>>>& graph);
	static std::map<Node, bool> invs(
			const std::map<Node, std::vector<std::pair<Node, bool>>>& graph);
	static void join(BoolGraph& graph1, const BoolGraph& graph2);
};
#endif
