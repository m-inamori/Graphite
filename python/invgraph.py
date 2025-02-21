from __future__ import annotations

# coding: utf-8
# invgraph.py

from typing import NewType, List, Tuple, Dict, Iterator

from graph import WeightedGraphBase, Node
from kruskal import Kruskal


#################### Type ####################

Edge = Tuple[Node, Node, float]
Value = Tuple[Node, float, bool]


#################### InvGraph ####################

class InvGraph(Dict[Node, List[Value]], WeightedGraphBase):
	def __init__(self) -> None:
		super().__init__(self)
		WeightedGraphBase.__init__(self)
	
	def __setitem__(self, key: Node, value: list[Value]) -> None:
		super().__setitem__(key, value)
	
	##### virtual methods for GraphBase #####
	def generate_nodes(self) -> Iterator[Node]:
		for v in self.keys():
			yield v
	
	def neighbors(self, v0: Node) -> list[Node]:
		value = self.get(v0, [])
		return [ v for v, d, b in value ]
	
	##### virtual methods for WeightedGraphBase #####
	def generate_weighted_edges(self) -> Iterator[Edge]:
		return ((u, v, w) for u, neis in self.items()
									for v, w, b in neis if u < v)
	
	##### non-virtual methods #####
	def minimum_spanning_tree(self) -> InvGraph:
		m: dict[tuple[Node, Node], Value] = { (u, v): (v, w, b)
												for u, neis in self.items()
												for v, w, b in neis }
		
		graph = InvGraph()
		for u in self.keys():
			graph[u] = []
		
		g: dict[Node, list[Node]] = Kruskal(self)
		for u, vs in g.items():
			for v in vs:
				_, w, b = m[(u, v)]
				graph[u].append((v, w, b))
		return graph
	
	def divide_into_connected(self) -> list[InvGraph]:
		nss: list[list[Node]] = self.divide_nodes_into_connected()
		graphs: list[InvGraph] = [ InvGraph() for _ in nss ]
		for i, ns in enumerate(nss):
			for n in ns:
				graphs[i][n] = self.__getitem__(n)
		return graphs
	
	def walk(self, v0: Node) -> list[Edge]:
		edges: list[Edge] = []
		visited: set[Node] = set()
		stack: list[Node] = [v0]
		visited.add(v0)
		while stack:
			v = stack.pop()
			vs = self.__getitem__(v)
			for v1, _, inv in vs:
				edges.append((v, v1, inv))
				if v1 in visited:
					continue
				stack.append(v1)
				visited.add(v1)
		return edges

__all__ = ['InvGraph']
