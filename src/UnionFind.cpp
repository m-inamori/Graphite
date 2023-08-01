#include "../include/UnionFind.h"

using namespace std;

UnionFind::UnionFind(const vector<Node>& nodes) {
	for(auto p = nodes.begin(); p != nodes.end(); ++p) {
		parents.insert(pair<size_t,size_t>(*p, *p));
		heights.insert(pair<size_t,int>(*p, 1));
	}
}

void UnionFind::join(const Node v1, const Node v2) {
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

UnionFind::Node UnionFind::root(Node v0) const {
	Node	v = v0;
	while(true) {
		auto	p = this->parents.find(v);
		if(p->second == v)
			return p->second;
		v = p->second;
	}
}

map<UnionFind::Node, size_t> UnionFind::divide_into_trees() const {
	map<Node, size_t>	groups;
	size_t	group_index = 0;
	for(auto p = parents.begin(); p != parents.end(); ++p) {
		const Node	v0 = p->first;
		auto	q = groups.find(v0);
		if(q != groups.end())	// defined
			continue;
		
		// find root or defined node
		vector<Node>	path(1, v0);
		while(true) {
			const Node	v = path.back();
			auto	r = groups.find(v);
			if(r != groups.end()) {	// defined
				break;
			}
			else {
				auto	p1 = parents.find(v);
				if(p1->second == v) {	// root
					groups.insert(make_pair(v, group_index));
					++group_index;
					break;
				}
				else {
					path.push_back(p1->second);
				}
			}
		}
		
		for(size_t j = 0; j < path.size() - 1; ++j)
			groups[path[j]] = groups[path.back()];
	}
	return groups;
}
