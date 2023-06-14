# coding: utf-8
# kruskal.py
# Kruskal { node: [(node, weight)] } -> { node: [(node, weight)] }

from __future__ import annotations
from itertools import *
from collections import defaultdict
from typing import Dict, List, Tuple, Generator, TypeVar


#################### UnionFind ####################

class UnionFind:
	def __init__(self, nodes: list[int]):
		self.parents: Dict[int, int] = { v: v for v in nodes }
		self.heights: Dict[int, int] = { v: 1 for v in nodes }
	
	def join(self, v1: int, v2: int):
		r1 = self.root(v1)
		r2 = self.root(v2)
		h1 = self.heights[r1]
		h2 = self.heights[r2]
		if h1 <= h2:
			self.parents[r1] = r2
			self.heights[r2] = max(self.heights[r2], self.heights[r1]+1)
		else:
			self.parents[r2] = r1
			self.heights[r1] = max(self.heights[r1], self.heights[r2]+1)
	
	def root(self, v: int) -> int:
		while self.parents[v] != v:
			v = self.parents[v]
		return v
	
	def __str__(self):
		def f(v: int) -> str:
			s = str(v)
			for child in children[v]:
				for line in f(child).split('\n'):
					s += '\n  ' + line
			return s
		
		children: Dict[int,list[int]] = { v: [] for v in self.parents.keys() }
		for v, parent in self.parents.items():
			if parent != v:
				children[parent].append(v)
		
		return '\n'.join(f(v) for v, parent in self.parents.items()
														if v == parent)


#################### UnionFindSerial ####################

class UnionFindSerial:
	def __init__(self, N: int):
		self.parents: list[int] = list(range(N))
		self.heights: list[int] = [1] * (N)
	
	def join(self, v1: int, v2: int):
		r1 = self.root(v1)
		r2 = self.root(v2)
		h1 = self.heights[r1]
		h2 = self.heights[r2]
		if h1 <= h2:
			self.parents[r1] = r2
			self.heights[r2] = max(self.heights[r2], self.heights[r1]+1)
		else:
			self.parents[r2] = r1
			self.heights[r1] = max(self.heights[r1], self.heights[r2]+1)
	
	def root(self, v: int) -> int:
		while self.parents[v] != v:
			v = self.parents[v]
		return v
	
	def __str__(self):
		def f(v: int) -> str:
			s = str(v)
			for child in children[v]:
				for line in f(child).split('\n'):
					s += '\n  ' + line
			return s
		
		children: Dict[int, list[int]] = { v: [] for v in self.parents.keys() }
		for v, parent in self.parents.items():
			if parent != v:
				children[parent].append(v)
		
		return '\n'.join(f(v) for v, parent in self.parents.items()
														if v == parent)


#################### process ####################

T = TypeVar('T', int, float)	# weight of edge of graph

def generate_edges(graph: Dict[int, list[tuple[int, T]]]
							) -> Generator[tuple[int, int, T], None, None]:
	for v1 in graph.keys():
		for v2, w in graph[v1]:
			if v1 < v2:
				yield (v1, v2, w)

def Kruskal(graph: Dict[int, list[tuple[int, T]]]
								) -> Dict[int, list[tuple[int, T]]]:
	nodes = list(graph.keys())
	tree = UnionFind(nodes)
	edges = sorted(generate_edges(graph), key=lambda v: v[2])
	
	new_graph = defaultdict(list)
	counter = 0
	for v1, v2, w in edges:
		if tree.root(v1) != tree.root(v2):
			tree.join(v1, v2)
			new_graph[v1].append((v2, w))
			new_graph[v2].append((v1, w))
			counter += 1
			if counter == len(nodes) - 1:
				break
	return new_graph

Graph = Dict[int, List[Tuple[int, T]]]

def Kruskal_fast(graph: Graph) -> Graph:
	dic: Dict[int, int] = dict(zip(graph.keys(), count()))
	N = len(dic)
	tree = UnionFindSerial(N)
	edges = sorted(generate_edges(graph), key=lambda v: v[2])
	
	new_graph: Graph = { v: [] for v in graph.keys() }
	counter = 0
	for v1, v2, w in edges:
		n1 = dic[v1]
		n2 = dic[v2]
		if tree.root(n1) != tree.root(n2):
			tree.join(n1, n2)
			new_graph[v1].append((v2, w))
			new_graph[v2].append((v1, w))
			counter += 1
			if counter == N - 1:
				break
	return new_graph


#################### main ####################

__all__ = ['Kruskal', 'Kruskal_fast']
