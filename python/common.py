from __future__ import annotations

# coding: utf-8
# common.py

from collections.abc import Iterable
from collections import defaultdict
import csv
import sys

from typing import TypeVar, Any, Generator, IO, Optional, Iterator, Dict


#################### I/O ####################

def read_csv(path: str, delim: str=',') -> Generator[list[str], None, None]:
	with open(path, 'r') as f:
		reader = csv.reader(f, delimiter=delim)
		for v in reader:
			yield v

def read_tsv(path: str) -> Generator[list[str], None, None]:
	for v in read_csv(path, '\t'):
		yield v

def write_csv(v: list[Any], out: IO):
	print(','.join(map(str_with_NA, v)), file=out)

def write_tsv(v: list[Any], out: IO):
	print('\t'.join(map(str_with_NA, v)), file=out)

def int_with_NA(s: str) -> Optional[int]:
	try:
		return int(s)
	except ValueError:
		return None

def str_with_NA(e) -> str:
	if e is None:
		return 'NA'
	else:
		return str(e)


#################### list ####################

def is_all_same(xs) -> bool:
	if all(x is None for x in xs):
		return True
	
	x0 = next(x for x in xs if x is not None)
	return all(x0 == x for x in xs)

def unique_list(*args):
	v = []
	for a in args:
		if a not in v:
			v.append(a)
	return v

T = TypeVar('T', bound=Iterable)

def flatten(v: T) -> list[Any]:
	w = []
	for e in v:
		if isinstance(e, Iterable):
			w.extend(flatten(e))
		else:
			w.append(e)
	return w


#################### graph ####################

def divide_graph_into_connected(g):
	def walk(v0):
		vs = [v0]
		stk = [v0]
		visited.add(v0)
		while stk:
			v1 = stk.pop()
			for e in g[v1]:
				v = e[0]
				if v in visited:
					continue
				vs.append(v)
				stk.append(v)
				visited.add(v)
		
		return vs
	
	if len(g) == 1:
		return [g]
	
	gs = []
	visited = set()
	for v in g.keys():
		if v in visited: continue
		vs = walk(v)
		gs.append(dict((w, g[w]) for w in vs if w in g))
	return gs

def nearest_vertecies(graph1, graph2):
	vs1 = list(graph1.keys())
	vs2 = list(graph2.keys())
	if vs1[-1] < vs2[0]:
		return (vs1[-1], vs2[0])
	elif vs2[-1] < vs1[0]:
		return (vs1[0], vs2[-1])
	
	d, v1, v2 = min((abs(v1 - v2), v1, v2) for v1, v2 in product(vs1, vs2))
	return (v1, v2)


#################### others ####################

K = TypeVar('K')
V = TypeVar('V')
def classify(iter: Iterator[tuple[K, V]]) -> Dict[K, list[V]]:
	dic: Dict[K, list[V]] = defaultdict(list)
	for k, v in iter:
		dic[k].append(v)
	return dic

