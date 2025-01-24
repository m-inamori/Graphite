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

class VCFHMM:
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	MIN_PROB = -1e300
	
	def __init__(self, records: list[VCFFamilyRecord],
								ref_haps: list[list[int]], map_: Map):
		self.ref_haps	= ref_haps
		self.map		= map_
		self.E = self.calc_E()
		self.Epc = self.calc_Epc()
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
		# 親は直接の後代でないので、乗り換え確率を高くする
		K = 5						# 親の遷移確率の子に対する倍率
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(records, records[1:]) ]
	
	def calc_E(self):
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E1 = [ [ 1.0-w*2, w/2,     w/2,     w ],
			   [ w/2,     1.0-w*2, w/2,     w ],
			   [ w,       1.0-3*w, w,       w ],
			   [ w/2,     w/2,     1.0-w*2, w ] ]
		
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	# 両親と後代のnon-phased genotypeからの後代の排出確率
	# あとで変える
	def calc_Epc(self):
		w = 0.01		# probability of wrong
		E1 = [ [ [ 1.0-w*3/4, w/2,       w/4,     w ],
				 [ 0.5-w*3/4, 0.5-w*3/4, w/2,     w ],
				 [ w/2,       1.0-w*2,   w/2,     w ],
				 [ 0.5-w*3/4, 0.5-w*3/4, w/2,     w ] ],
			   [ [ 0.5-w*3/4, 0.5-w*3/4, w/2,     w ],
			     [ 1/4-w/3,   0.5-w/3,   1/4-w/3, w ],
			     [ w,       0.5-w,   0.5-w,   w ],
			     [ 1/4-w/3, 0.5-w/3, 1/4-w/3, w ] ],
			   [ [ w,       1.0-3*w, w,       w ],
			     [ w,       0.5-w,   0.5-w,   w ],
			     [ w,       w,       1.0-3*w, w ],
			     [ w,       0.5-w,   0.5-w,   w ] ] ]
		E = [ [ [ log(p) for p in v ] for v in w ] for w in E1 ]
		return E
	
	def parent_transition_probability(self, i: int, prev_hp: int, hp) -> float:
		cp = self.Cp[i-1]	# 親の遷移確率
		hp2, hp1 = divmod(h, NH)
		prev_hp2, prev_hp1 = divmod(prev_h, NH)
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cp if prev_hp1 != hp1 else 1.0 - cp) +
				log(cp if prev_hp2 != hp2 else 1.0 - cp))
	
	def progeny_transition_probability(self, i: int, prev_hc: int, hc) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		hc1 = hc & 1
		hc2 = hc >> 1
		prev_hc1 = prev_hc & 1
		prev_hc2 = prev_hc >> 1
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
				log(cc if prev_hc2 != hc2 else 1.0 - cc))
	
	def trace_back(dps: list[DP]) -> list[int]:
		hs = [ 0 ] * M
		h, (_, prev_h) = max(enumerate(dp[-1]), key=lambda v: v[1])
		hs[-1] = h
		hs[-2] = prev_h
		for i in range(M-2, 0, -1):
			_, prev_h = dp[i][prev_h]
			hs[i-1] = prev_h
		return hs
	
	def compute_phased_gt_by_refhaps(hp: int, i: int) -> int:
		NH = len(self.ref_haps)
		hp2, hp1 = divmod(hp, NH)
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
