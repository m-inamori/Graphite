from __future__ import annotations

# coding: utf-8
# graph.py

from abc import ABC, abstractmethod
from typing import TypeVar, NewType, Tuple, Iterator, Union

from UnionFind import *


#################### Type ####################

Edge = Tuple[Node, Node]
Weight = Union[int, float]
WeightedEdge = Tuple[Node, Node, Weight]


#################### GraphBase ####################

class GraphBase(ABC):
	def __init__(self):
		pass
	
	@abstractmethod
	def generate_nodes(self) -> Iterator[Node]:
		pass
	
	@abstractmethod
	def neighbors(self, v0: Node) -> list[Node]:
		pass
	
	def generate_edges(self) -> Iterator[Edge]:
		for u in self.generate_nodes():
			for v in self.neighbors(u):
				yield (u, v)
	
	def divide_nodes_into_connected(self) -> list[list[Node]]:
		nodes = list(self.generate_nodes())
		uf = UnionFind(nodes)
		for u, v in self.generate_edges():
			uf.join(u, v)
		
		groups: dict[Node, int] = uf.divide_into_trees()
		num_groups = max(groups.values()) + 1
		vs: list[list[Node]] = [ [] for _ in range(num_groups) ]
		for v, group_index in groups.items():
			vs[group_index].append(v)
		return vs


#################### WeightedGraphBase ####################

class WeightedGraphBase(GraphBase):
	def __init__(self):
		pass
	
	@abstractmethod
	def generate_weighted_edges(self) -> Iterator[WeightedEdge]:
		pass

__all__ = ['GraphBase', 'WeightedGraphBase', 'Node']
