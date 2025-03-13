#ifndef __UNIONFIND
#define __UNIONFIND

#include <vector>
#include <map>


class UnionFind {
public:
	using Node = std::size_t;
	
private:
	std::map<Node, Node>	parents;
	std::map<Node, int>		heights;
	
public:
	explicit UnionFind(const std::vector<Node>& nodes);
	
	void join(Node v1, Node v2);
	Node root(Node v0) const;
	std::map<Node, std::size_t> divide_into_trees() const;
};
#endif
