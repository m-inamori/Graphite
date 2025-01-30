from __future__ import annotations

# coding: utf-8
# VCFHMM.py
# 補完するためにHMMを利用する関数を提供するVCF

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple

from VCFFamily import *
from Map import *
from Genotype import Genotype


#################### VCFHMM ####################

class VCFHMM(VCFMeasurable):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord], map_: Map):
		VCFMeasurable.__init__(self, map_)
		self.map		= map_
		self.E = self.calc_E()
	
	def calc_E(self):
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E1 = [ [ 1.0-w*2, w/2,     w/2,     w ],
			   [ w/2,     1.0-w*2, w/2,     w ],
			   [ w/2,     1.0-w*2, w/2,     w ],
			   [ w/2,     w/2,     1.0-w*2, w ] ]
		
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	def dist(self, r1: VCFRecord, r2: VCFRecord) -> float:
		dist = (self.cM(r2.pos()) - self.cM(r1.pos())) / 100
		if dist > 0.0:
			return dist
		else:
			return (r2.pos() - r1.pos()) * 1e-6
	
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
