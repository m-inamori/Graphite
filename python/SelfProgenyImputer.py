from __future__ import annotations

# coding: utf-8
# SelfProgenyImputer.py
# 親が自殖で親がimputedなとき後代をimputeする

from math import log
from typing import List, Tuple

from VCF import *
from VCFHMM import *
from Genotype import Genotype


#################### SelfProgenyImputer ####################

MIN_PROB = -1e300

class SelfProgenyImputer(VCFHMM[VCFRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFRecord],
							iprog: int, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFRecord] = records
		self.ic = iprog		# 何番目の後代か
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def num_progenies(self) -> int:
		return len(self.records[0].v) - 10
	
	def progeny_genotype(self, h: int, parent_gt: int) -> int:
		hc2, hc1 = divmod(h, 2)
		return ((parent_gt >> hc1) & 1) | (((parent_gt >> hc2) & 1) << 1)
	
	# i: record index
	def emission_probability(self, h: int, parent_gt: int, oc: int) -> float:
		gt_prog = self.progeny_genotype(h, parent_gt)
		return self.E[gt_prog][oc]
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		dp = [ [ (MIN_PROB, 0) ] * 4 for _ in range(M) ]
		record = self.records[0]
		parent_gt = Genotype.phased_gt_to_int(self.records[0].v[9])
		# observed progeny
		oc = Genotype.phased_gt_to_int(record.v[self.ic+10])
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
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		record = self.records[i]
		cc = self.Cc[i-1]	# 遷移確率
		parent_gt = Genotype.phased_gt_to_int(self.records[i].v[9])
		# observed progs
		oc = Genotype.gt_to_int(record.v[self.ic+10])
		
		for h in range(4):		# hidden state
			E_all = self.emission_probability(h, parent_gt, oc)
			
			for prev_h in range(4):
				T_all = self.transition_probability(h, prev_h, cc)
				
				prob = dp[i-1][prev_h][0] + (T_all + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		N = self.num_progenies()
		for i in range(M):
			record = self.records[i]
			parent_gt = Genotype.phased_gt_to_int(record.v[9])
			prog_gt = self.progeny_genotype(hs[i], parent_gt)
			record.set_GT(self.ic + 1, Genotype.int_to_phased_gt(prog_gt))
	
	def impute(self) -> None:
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
