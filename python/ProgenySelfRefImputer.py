from __future__ import annotations

# coding: utf-8
# ProgenyRefImputer.py
# 親がphasedのときに後代をimputeする
# リファレンスがあるとき用
# phasedの遺伝子型を考慮する

from math import log
from typing import List, Tuple

from GenoRecord import GenoRecord
from VCFHMMSelfRef import *
from Genotype import Genotype


#################### ProgenySelfRefImputer ####################

MIN_PROB = -1e300

class ProgenySelfRefImputer(VCFHMMSelfRef[GenoRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[GenoRecord],
									map_: Map, w: float) -> None:
		VCFHMMSelfRef.__init__(self, records, map_)
		self.records: list[GenoRecord] = records
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 2
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc: int, parent_gt: int) -> int:
		hc1 = hc & 1
		hc2 = hc >> 1
		return ((parent_gt >> hc1) & 1) | (((parent_gt >> hc2) & 1) << 1)
	
	# i: record index, j: progeny index
	def emission_probability(self, i: int, j: int, h: int,
											parent_gt: int) -> float:
		record = self.records[i]
		gt = record.geno[j+1]
		oc = 4 if gt == Genotype.NA else gt & 3
		phased_gt = self.gt_by_haplotypes(h, parent_gt)
		return self.E[phased_gt][oc]
	
	def progeny_transition_probability(self, i: int,
											prev_hc: int, hc: int) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		hc1 = hc & 1
		hc2 = hc >> 1
		prev_hc1 = prev_hc & 1
		prev_hc2 = prev_hc >> 1
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
				log(cc if prev_hc2 != hc2 else 1.0 - cc))
	
	def initialize_dp(self, j: int, M: int) -> list[DP]:
		dp = [ [ (MIN_PROB, 0) ] * 4 for _ in range(M) ]
		record = self.records[0]
		parent_gt = record.geno[0]
		for h in range(4):		# hidden state
			E = self.emission_probability(0, j, h, parent_gt)
			dp[0][h] = (E, h)
		return dp
	
	def update_dp(self, i: int, j: int, dp: list[DP]) -> None:
		record = self.records[i]
		parent_gt = record.geno[0] & 3
		
		for hc in range(4):		# hidden state
			E_all = self.emission_probability(i, j, hc, parent_gt)
			
			cc = self.Cc[i-1]	# transition probability
			hc1 = hc & 1
			hc2 = hc >> 1
			
			for prev_hc in range(4):
				prev_hc1 = prev_hc & 1
				prev_hc2 = prev_hc >> 1
				# 遷移確率 0なら(1-c)、1ならcを掛ける
				Tc = (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
					  log(cc if prev_hc2 != hc2 else 1.0 - cc))
				
				prob = dp[i-1][prev_hc][0] + (Tc + E_all)
				dp[i][hc] = max(dp[i][hc], (prob, prev_hc))
	
	def update_genotypes(self, j: int, hs: list[int]) -> None:
		M = len(self.records)
		for i in range(M):
			record = self.records[i]
			parent_gt = record.geno[0]
			gtc_int = self.gt_by_haplotypes(hs[i], parent_gt)
			record.geno[j+1] = gtc_int | 4
	
	# j: index of progeny
	def impute(self, j: int) -> None:
		M = len(self.records)	# マーカー数
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = self.initialize_dp(j, M)
		for i in range(1, M):
			self.update_dp(i, j, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(j, hs)
