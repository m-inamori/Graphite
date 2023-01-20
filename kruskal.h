#ifndef __KRUSKAL
#define __KRUSKAL

#include <string>
#include <vector>
#include <map>


namespace Kruskal {
	
	using GRAPH = std::map<std::size_t,
							std::vector<std::pair<std::size_t,double>>>;
	
	
	//////////////////// UnionFind ////////////////////
	
	class UnionFind {
		std::map<std::size_t,std::size_t>	parents;
		std::map<std::size_t,int>	heights;
		
	public:
		UnionFind(const std::vector<std::size_t>& nodes);
		
		void join(const std::size_t& v1, const std::size_t& v2);
		const std::size_t& root(const std::size_t& v0) const;
	};
	
	//////////////////// Kruskal ////////////////////
	
	GRAPH Kruskal(const GRAPH& graph);
}
#endif
