from __future__ import annotations

# coding: utf-8
# inverse_graph.py
# Nodeを反転すべきかどうかを決めるグラフ

from itertools import count, product
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from option import *
from common import flatten, divide_graph_into_connected

def is_consistent(graph: Dict[int, list[tuple[int, bool]]]) -> bool:
	visited: Dict[int, bool] = { }
	v0 = next(v for v in graph.keys())
	stk = [(v0, False)]		# v0は反転しないとする
	while stk:
		v, b = stk.pop()
		for v1, b1_edge in graph[v]:
			b1 = b ^ b1_edge
			if v1 in visited:
				b1_prev = visited[v1]
				if b1 != b1_prev:
					return False
			else:
				visited[v1] = b1
				stk.append((v1, b1))
	return True
	
def invs(graph: Dict[int, list[tuple[int, bool]]]) -> Dict[int, bool]:
	visited: Dict[int, bool] = { }
	v0 = next(v for v in graph.keys())
	stk = [(v0, False)]		# v0は反転しないとする
	while stk:
		v, b = stk.pop()
		for v1, b1_edge in graph[v]:
			b1 = b ^ b1_edge
			if v1 not in visited:
				visited[v1] = b1
				stk.append((v1, b1))
	return visited


#################### InverseGraph ####################

class InverseGraph:
	def __init__(self, g: Dict[int, list[tuple[int, int, int]]]):
		self.g: Dict[int, list[tuple[int, int, int]]] = g
	
	def __len__(self):
		return len(self.g)
	
	def optimize_inversions(self) -> tuple[bool, ...]:
		subgraphs = self.divide_graph_into_connected()
		bs = [False] * len(self)
		for subg in subgraphs:
			dic_bs = subg.optimize_inversions_connected()
			for i, b in dic_bs.items():
				bs[i] = b
		return tuple(bs)
	
	def divide_graph_into_connected(self) -> list[InverseGraph]:
		subgraphs = divide_graph_into_connected(self.g)
		return [ InverseGraph(g) for g in subgraphs ]
	
	# selfは連結成分である前提
	def optimize_inversions_connected(self) -> Dict[int, bool]:
		if len(self) <= 20:
			return self.search_all()
		
		# trueとfalseの偏りが大きいものからedgeを確定させていく
		# グラフが全て繋がるまでに矛盾が無ければOK
		dic_bs = self.connect_biased_edges()
		if dic_bs:
			return dic_bs
		
		# それでもダメならランダム
		# 本当はSAを使いたい
		return self.search_randomly()
	
	# しらみつぶし
	def search_all(self) -> Dict[int, bool]:
		if len(self) == 1:
			return { list(self.g.keys())[0]: False }
		
		bs = min(product((False, True), repeat=len(self)-1),
										key=lambda bs: self.score(bs))
		vs = sorted(self.g.keys())
		return dict(zip(vs, (False,) + bs))
	
	def search_randomly(self) -> Dict[int, bool]:
		bs = min((tuple(random.choices((False, True), k=len(self)-1))
												for _ in range(2**19)),
												key=lambda bs: self.score(bs))
		vs = sorted(self.g.keys())
		return dict(zip(vs, (False,) + bs))
	
	def generate_edges(self) -> Iterator[tuple[int, int, int, int]]:
		for i, v in self.g.items():
			for j, n1, n2 in v:
				if i > j:
					continue
				yield (i, j, n1, n2)
	
	def score(self, bs: tuple[bool, ...]) -> int:
		dic: Dict[int, int] = dict(zip(sorted(self.g.keys()), count()))
		num_wrong: int = 0
		for v1, v2, n1, n2 in self.generate_edges():
			i = dic[v1]
			j = dic[v2]
			b1 = False if i == 0 else bs[i-1]
			b2 = bs[j-1]
			if b1 == b2:
				num_wrong += n2
			else:
				num_wrong += n1
		return num_wrong
	
	def connect_biased_edges(self) -> Dict[int, bool]:
		def ratio(v: tuple[int, int, int, int]) -> float:
			i, j, n1, n2 = v
			return min((n1 + 1) / (n1 + n2 + 2), (n2 + 1) / (n1 + n2 + 2))
		
		def join(graph1: Dict[int, list[tuple[int, bool]]],
							graph2: Dict[int, list[tuple[int, bool]]]):
			for v1, vs in graph2.items():
				graph1[v1] = vs
		
		sorted_edges = sorted(((ratio(e), e) for e in self.generate_edges()))
		subgraphs: list[Dict[int, list[tuple[int, bool]]]] = \
									[ { v: [] } for v in self.g.keys() ]
		for r, (i, j, n1, n2) in sorted_edges:
			k1 = next(k for k, g in enumerate(subgraphs) if i in g.keys())
			k2 = next(k for k, g in enumerate(subgraphs) if j in g.keys())
			if k1 == k2:
				g1 = subgraphs[k1]
				g1[i].append((j, n1 < n2))
				g1[j].append((i, n1 < n2))
				if not is_consistent(g1):
					return { }
			else:
				subg1 = subgraphs[k1]
				subg2 = subgraphs.pop(k2)
				join(subg1, subg2)
				subg1[i].append((j, n1 < n2))
				subg1[j].append((i, n1 < n2))
				if len(subgraphs) == 1:
					return invs(subgraphs[0])
		else:
			return { }	# ここには来ないはず
	
	@staticmethod
	def convert(graph: list[list[tuple[int, int, int]]]) -> InverseGraph:
		return InverseGraph(dict(enumerate(graph)))


__all__ = ['InverseGraph']
