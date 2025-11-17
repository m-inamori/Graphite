# coding: utf-8
# kruskal.py
# Kruskal { node: [(node, weight)] } -> { node: [(node, weight)] }

from __future__ import annotations
from itertools import *
from collections import defaultdict
from typing import Iterator, TypeVar, Union

from graph import WeightedGraphBase, Node
from UnionFind import *


#################### process ####################

Weight = Union[int, float]	# weight of edge of graph

def Kruskal(graph: WeightedGraphBase) -> dict[Node, list[Node]]:
	nodes = list(graph.generate_nodes())
	tree = UnionFind(nodes)
	edges = sorted(graph.generate_weighted_edges(),
									key=lambda v: (v[2], v[0], v[1]))
	
	new_graph: dict[Node, list[Node]] = defaultdict(list)
	counter = 0
	for v1, v2, w in edges:
		if tree.root(v1) != tree.root(v2):
			tree.join(v1, v2)
			new_graph[v1].append(v2)
			new_graph[v2].append(v1)
			counter += 1
			if counter == len(nodes) - 1:
				break
	return new_graph


#################### main ####################

__all__ = ['Kruskal']
