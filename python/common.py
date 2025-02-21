# coding: utf-8
# common.py

from __future__ import annotations
from itertools import product
from collections import defaultdict
import csv
import sys

from typing import TypeVar, Any, Optional, Iterator, Sequence, TextIO


#################### I/O ####################

def read_csv(path: str, delim: str = ',') -> Iterator[list[str]]:
	with open(path, 'r') as f:
		reader = csv.reader(f, delimiter=delim)
		for v in reader:
			yield v

def read_tsv(path: str) -> Iterator[list[str]]:
	for v in read_csv(path, '\t'):
		yield v

# separeted by white space
def read_wsv(path: str) -> Iterator[list[str]]:
	with open(path, 'r') as f:
		for line in f:
			yield line.split()

def write_csv(v: list[Any], out: TextIO) -> None:
	print(','.join(map(str_with_NA, v)), file=out)

def write_tsv(v: list[Any], out: TextIO) -> None:
	print('\t'.join(map(str_with_NA, v)), file=out)

def int_with_NA(s: int) -> Optional[int]:
	try:
		return int(s)
	except ValueError:
		return None

def str_with_NA(e: Any) -> str:
	if e is None:
		return 'NA'
	else:
		return str(e)

T = TypeVar('T')

def is_all_same(xs: Sequence[T]) -> bool:
	if all(x is None for x in xs):
		return True
	
	x0 = next(x for x in xs if x is not None)
	return all(x0 == x for x in xs)

def unique_list(*args: T) -> list[T]:
	v = []
	for a in args:
		if a not in v:
			v.append(a)
	return v

K = TypeVar('K')
V = TypeVar('V')

def classify(iterable: Iterator[tuple[K, V]]) -> dict[K, list[V]]:
	dic = defaultdict(list)
	for k, v in iterable:
		dic[k].append(v)
	return dic
