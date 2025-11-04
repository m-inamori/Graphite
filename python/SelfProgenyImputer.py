from __future__ import annotations

# coding: utf-8
# SelfProgenyImputer.py
# 親が自殖で親がimputedなとき後代をimputeする

from math import log
from typing import List, Tuple

from GenoRecord import GenoRecord
from VCFHMM import *
from Genotype import Genotype


#################### SelfProgenyImputer ####################

MIN_PROB = -1e300

class SelfProgenyImputer(VCFHMM[GenoRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[GenoRecord], map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[GenoRecord] = records
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 1
	
	def progeny_genotype(self, h: int, parent_gt: int) -> int:
		hc2, hc1 = divmod(h, 2)
		return ((parent_gt >> hc1) & 1) | (((parent_gt >> hc2) & 1) << 1)
	
	# i: record index
	def emission_probability(self, h: int, parent_gt: int, oc: int) -> float:
		gt_prog = self.progeny_genotype(h, parent_gt)
		return self.E[gt_prog][oc&3]
	
	def initialize_dp(self, iprog: int) -> list[DP]:
		M = len(self.records)	# マーカー数
		dp = [ [ (MIN_PROB, 0) ] * 4 for _ in range(M) ]
		record = self.records[0]
		parent_gt = self.records[0].geno[0]
		# observed progeny
		oc = record.geno[iprog+1]
		for h in range(4):		# hidden state
			E_all = self.emission_probability(h, parent_gt, oc)
			dp[0][h] = (E_all, h)
		return dp
	
	def transition_probability(self, h: int, prev_h: int, cc: float) -> float:
		h2, h1 = divmod(h, 2)
		prev_h2, prev_h1 = divmod(prev_h, 2)
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if h1 != prev_h1 else 1.0 - cc) +
				log(cc if h2 != prev_h2 else 1.0 - cc))
	
	def update_dp(self, i: int, iprog: int, dp: list[DP]) -> None:
		record = self.records[i]
		cc = self.Cc[i-1]	# 遷移確率
		parent_gt = self.records[i].geno[0]
		# observed progs
		oc = record.unphased(iprog+1)
		
		for h in range(4):		# hidden state
			E_all = self.emission_probability(h, parent_gt, oc)
			
			for prev_h in range(4):
				T_all = self.transition_probability(h, prev_h, cc)
				
				prob = dp[i-1][prev_h][0] + (T_all + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int], iprog: int) -> None:
		M = len(self.records)
		N = self.num_progenies()
		for i in range(M):
			record = self.records[i]
			parent_gt = record.geno[0]
			prog_gt = self.progeny_genotype(hs[i], parent_gt)
			record.geno[iprog+1] = prog_gt | 4
	
	def impute(self, iprog: int) -> None:
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		M = len(self.records)
		dp = self.initialize_dp(iprog)
		for i in range(1, M):
			self.update_dp(i, iprog, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs, iprog)
