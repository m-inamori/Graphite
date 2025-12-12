from __future__ import annotations

# coding: utf-8
# ProgenyImputerByOneParent.py
# 片親がphasedでもう片親がunknownのときに後代をimputeする

from math import log
from typing import List, Tuple

from VCFFamily import VCFFamilyRecord
from VCFHMM import *
from Genotype import Genotype


#################### ProgenyImputerByOneParent ####################

MIN_PROB = -1e300

class ProgenyImputerByOneParent(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
						ref_haps: list[list[int]],
						is_mat_imputed: bool, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_imputed = is_mat_imputed
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
		K = 5						# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(records, records[1:]) ]
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 2
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def parent_alleles(self, h: int, i: int) -> tuple[int, int]:
		hc1 = h & 1
		hc2 = h >> 1
		a2 = self.ref_haps[hc2][i]
		if self.is_mat_imputed:
			a1 = (self.records[i].geno[0] >> hc1) & 1
			return (a1, a2)
		else:
			a1 = (self.records[i].geno[1] >> hc1) & 1
			return (a2, a1)
	
	# i: record index, j: progeny index
	def emission_probability(self, i: int, j: int, h: int,
											parent_gt: int) -> float:
		oc = self.records[i].unphased(j+2)
		mat_a, pat_a = self.parent_alleles(h, i)
		phased_gt = Genotype.from_alleles(mat_a, pat_a) & 3
		return self.E[phased_gt][oc]
	
	def transition_probability(self, i: int, prev_h: int, h: int) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		cp = self.Cp[i-1]
		hc1 = h & 1
		hc2 = h >> 1
		prev_hc1 = prev_h & 1
		prev_hc2 = prev_h >> 1
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
				log(cp if prev_hc2 != hc2 else 1.0 - cp))
	
	def initialize_dp(self, j: int, M: int) -> list[DP]:
		L = self.NH() * 2
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		parent_gt = record.mat_gt() if self.is_mat_imputed else record.pat_gt()
		for h in range(L):		# hidden state
			E = self.emission_probability(0, j, h, parent_gt)
			dp[0][h] = (E, h)
		return dp
	
	def update_dp(self, i: int, j: int, dp: list[DP]) -> None:
		L = self.NH() * 2
		record = self.records[i]
		parent_gt = record.mat_gt() if self.is_mat_imputed else record.pat_gt()
		
		for h in range(L):		# hidden state
			E = self.emission_probability(i, j, h, parent_gt)
			for prev_h in range(L):
				T = self.transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (T + E)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, j: int, hs: list[int]) -> None:
		M = len(self.records)
		for i in range(M):
			record = self.records[i]
			mat_a, pat_a = self.parent_alleles(hs[i], i)
			record.geno[j+2] = Genotype.from_alleles(mat_a, pat_a)
	
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
