from __future__ import annotations

# coding: utf-8
# VCFHMMRef.py
# 補完するためにHMMを利用する関数を提供するVCF
# リファレンスで補完するためのサブクラス

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple, TypeVar, Generic, Sequence

from GenoRecord import GenoRecord
from Map import *
from Genotype import Genotype


#################### VCFHMMRef ####################

R = TypeVar('R', bound=GenoRecord)

class VCFHMMRef(Generic[R], VCFMeasurable):
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
		E1 = [ [ 1.0-w*2, w/3,         w/3,         w/3,     w ],
			   [ w/3,     0.9-5.3/3*w, 0.1+0.1*w,   w/3,     w ],
			   [ w/3,     0.1+0.1*w,   0.9-5.3/3*w, w/3,     w ],
			   [ w/3,     w/3,         w/3,         1.0-w*2, w ] ]
		
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
