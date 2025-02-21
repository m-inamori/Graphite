from __future__ import annotations

# coding: utf-8
# ProgenyImputer.py
# 親がphasedのときに後代をimputeする

from math import log
from typing import List, Tuple

from VCFHMM import *
from Genotype import Genotype


#################### ProgenyImputer ####################

MIN_PROB = -1e300

class ProgenyImputer(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
									map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
	
	def num_progenies(self) -> int:
		return len(self.records[0].v) - 11
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc: int, mat_gt: int, pat_gt: int) -> int:
		hc1 = hc & 1
		hc2 = hc >> 1
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	# i: record index, j: progeny index
	def emission_probability(self, i: int, j: int, h: int,
									mat_gt: int, pat_gt: int) -> float:
		record = self.records[i]
		oc = Genotype.gt_to_int(record.v[j+11])
		phased_gt = self.gt_by_haplotypes(h, mat_gt, pat_gt)
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
		mat_gt = Genotype.phased_gt_to_int(record.v[9])
		pat_gt = Genotype.phased_gt_to_int(record.v[10])
		for h in range(4):		# hidden state
			E = self.emission_probability(0, j, h, mat_gt, pat_gt)
			dp[0][h] = (E, h)
		return dp
	
	def update_dp(self, i: int, j: int, dp: list[DP]) -> None:
		record = self.records[i]
		mat_gt = Genotype.phased_gt_to_int(record.v[9])
		pat_gt = Genotype.phased_gt_to_int(record.v[10])
		
		for hc in range(4):		# hidden state
			E_all = self.emission_probability(i, j, hc, mat_gt, pat_gt)
			
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
			mat_gt = Genotype.phased_gt_to_int(record.v[9])
			pat_gt = Genotype.phased_gt_to_int(record.v[10])
			gtc_int = self.gt_by_haplotypes(hs[i], mat_gt, pat_gt)
			record.v[j+11] = Genotype.int_to_phased_gt(gtc_int)
	
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
