from __future__ import annotations

# coding: utf-8
# ClassifyRecord.py

from functools import reduce
from collections import defaultdict

from VCFImpFamily import *
from VCFHomoHomo import *
from VCFHeteroHomo import *
from VCFHeteroHeteroLite import *
from pedigree import *
from TypeDeterminer import *
from common import *


#################### classify ####################

memo_tds: Dict[tuple[int, float], TypeDeterminer] = { }

def count_int_gts(gts: list[int]) -> tuple[int,int,int]:
	ns = [0, 0, 0]
	for gt in gts:
		if gt != -1:
			ns[gt] += 1
	return (ns[0], ns[1], ns[2])

def classify_record(record: VCFRecord, td: TypeDeterminer) -> tuple[str, int]:
	gts = record.get_int_gts()
	counter = count_int_gts(gts[2:])
	pairs = td.determine(counter)
	pair, wrong_type = classify_record_core(pairs, gts[0], gts[1])
	return (wrong_type, pair)

def classify_record_core(pairs: list[tuple[int, float]],
								mat_gt: int, pat_gt: int) -> tuple[int, str]:
	def is_matched(mat_gt: int, pat_gt: int, pair: int) -> bool:
		gt_pair = TypeDeterminer.int_gt_pair(pair)
		return (mat_gt, pat_gt) == gt_pair or (pat_gt, mat_gt) == gt_pair
	
	# pairsの中に確率が大きいものがあり、それ以外は小さければ、それだけにする
	if len(pairs) >= 2:
		pairs.sort(key=lambda v: v[1])	# 1-確率でソート
		ps = [ p for pair, p in pairs ]
		P = [ [ 1.0-ps[j] if j == i else ps[j]+0.01 for j in range(len(pairs)) ]
													for i in range(len(pairs)) ]
		qs = [ reduce(lambda x, y: x * y, v, 1.0) for v in P ]
		pair, p = pairs[0]
		if qs[0] / sum(qs) >= 0.99:
			pairs = pairs[:1]
	
	if not pairs:
		return (-1, 'Unspecified')
	elif len(pairs) == 1:
		pair, p = pairs[0]
		if pair in (0, 2, 5):	# same parent gt
			gt = pair >> 1
			if mat_gt == gt and pat_gt == gt:
				return (pair, 'Right')
			else:
				return (pair, 'Modifiable')
		else:
			# 0/0 x 0/1 -> 2, 0/0 x 1/1 -> 1, 0/1 x 1/1 -> 0
			avoiding_gt = (5 - pair) >> 1
			if mat_gt == pat_gt:
				return (pair, 'Unmodifiable')
			elif mat_gt == -1 and pat_gt not in (-1, avoiding_gt):
				return (pair, 'Modifiable')
			elif pat_gt == -1 and mat_gt not in (-1, avoiding_gt):
				return (pair, 'Modifiable')
			elif mat_gt != avoiding_gt and pat_gt != avoiding_gt:
				return (pair, 'Right')
			else:
				return (pair, 'Modifiable')
	else:
		pair, p = pairs[0]
		if is_matched(mat_gt, pat_gt, pair):
			return (pair, 'Right')
		if mat_gt == pat_gt:
			return (-1, 'Mix')
		
		# 最も優先順位が高いペアだけ片方だけ合っているなら修正可能
		bs = []
		for pair, p in pairs:
			gt_pair = TypeDeterminer.int_gt_pair(pair)
			bs.append(mat_gt in gt_pair or pat_gt in gt_pair)
		
		if bs[0] and all(not b for b in bs[1:]):
			return (pairs[0][0], 'Modifiable')
		else:
			return (-1, 'Mix')

def prepare(n: int, p: float):
	if (n, p) not in memo_tds:
		td = TypeDeterminer(n, p)
		memo_tds[(n, p)] = td

def get_typedeterminer(n: int, p: float) -> TypeDeterminer:
	if n in memo_tds:
		return memo_tds[(n, p)]
	
	td = TypeDeterminer(n, p)
	memo_tds[(n, p)] = td
	return td


#################### main ####################

__all__ = ['classify_record', 'get_typedeterminer', 'prepare']
