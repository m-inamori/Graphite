from __future__ import annotations

# coding: utf-8
# VCFHMMSelfRef.py
# 補完するためにHMMを利用する関数を提供するVCF
# リファレンスで補完するためのサブクラス

from itertools import product
from math import log
from typing import List, Tuple, TypeVar, Generic, Sequence

from GenoRecord import GenoRecord
from Map import *
from Genotype import Genotype


#################### VCFHMMSelfRef ####################

R = TypeVar('R', bound=GenoRecord)

class VCFHMMSelfRef(Generic[R], VCFMeasurable):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: Sequence[GenoRecord], map_: Map):
		VCFMeasurable.__init__(self, map_)
		self.map		= map_
		self.E = self.calc_E()
	
	def calc_E(self) -> list[list[float]]:
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3, N/A: 4
		# 0|1 <=> 1|0は0.1の確率で間違えるとする
		w = 0.01		# probability of wrong
		E1 = [ [0.0] * 5 for _ in range(4) ]
		for gt in range(4):
			orig_prob = [0.0] * 5
			for i, j in product(range(2), repeat=2):
				a1 = (gt >> i) & 1
				a2 = (gt >> j) & 1
				prog_gt = a1 | (a2 << 1)
				if prog_gt == 1:
					orig_prob[1] += 0.25 * 0.9
					orig_prob[2] += 0.25 * 0.1
				elif prog_gt == 2:
					orig_prob[1] += 0.25 * 0.1
					orig_prob[2] += 0.25 * 0.9
				else:
					orig_prob[gt] += 0.25
			# エラーを入れる
			for prog_gt in range(4):
				for error_gt in range(5):
					if error_gt == 4:
						E1[gt][error_gt] = w
					elif error_gt == prog_gt:
						E1[gt][error_gt] += orig_prob[gt] * (1.0 - w*2)
					else:
						E1[gt][error_gt] += orig_prob[gt] * (w/3)
		
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	def dist(self, r1: GenoRecord, r2: GenoRecord) -> float:
		dist = (self.cM(r2.pos) - self.cM(r1.pos)) / 100
		if dist > 0.0:
			return dist
		else:
			return (r2.pos - r1.pos) * 1e-6
	
	def trace_back(self, dp: list[DP]) -> list[int]:
		M = len(dp)
		hs: list[int] = [ 0 ] * M
		h, (_, prev_h) = max(enumerate(dp[-1]), key=lambda v: v[1])
		hs[-1] = h
		hs[-2] = prev_h
		for i in range(M-2, 0, -1):
			_, prev_h = dp[i][prev_h]
			hs[i-1] = prev_h
		return hs
