from __future__ import annotations

# coding: utf-8
# OrphanImputer.py
# 孤立したサンプルをimputeする

from math import log
from typing import List, Tuple

from GenoRecord import GenoRecord
from VCFHMM import *
from Genotype import Genotype


#################### OrphanImputer ####################

MIN_PROB = -1e300

class OrphanImputer(VCFHMM[GenoRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[GenoRecord],
						ref_haps: list[list[int]], map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[GenoRecord] = records
		self.ref_haps	= ref_haps
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		# crossover values
		# 親は直接の後代でないので、乗り換え確率を高くする
		K = 5						# 親の遷移確率の子に対する倍率
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def M(self) -> int:
		return len(self.ref_haps[0])
	
	def compute_phased_gt_by_refhaps(self, hp: int, i: int) -> int:
		hp2, hp1 = divmod(hp, self.NH())
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def emission_probability(self, i: int, h: int, orphan_gt: int) -> float:
		record = self.records[i]
		phased_gt = self.compute_phased_gt_by_refhaps(h, i)
		return self.E[phased_gt][orphan_gt]
	
	def orphan_transition_probability(self, i: int, prev_hp: int, hp: int) -> float:
		cp = self.Cp[i-1]	# 親の遷移確率
		hp2, hp1 = divmod(hp, self.NH())
		prev_hp2, prev_hp1 = divmod(prev_hp, self.NH())
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cp if prev_hp1 != hp1 else 1.0 - cp) +
				log(cp if prev_hp2 != hp2 else 1.0 - cp))
	
	def initialize_dp(self, io: int, L: int, M: int) -> list[DP]:
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		orphan_gt = record.unphased(io)
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h, orphan_gt)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		L = self.NH()**2
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			hp2, hp1 = divmod(h, self.NH())
			# 両側乗り換えることはないとする
			# non-phasedの親のあり得る前の状態
			for h1 in range(self.NH()):
				prev_h1 = h1 + hp2 * self.NH()
				prev_h_table[h].append(prev_h1)
			for h2 in range(self.NH()):
				if h2 == hp2:	# for duplication
					continue
				prev_h1 = hp1 + h2 * self.NH()
				prev_h_table[h].append(prev_h1)
		
		return prev_h_table
	
	def update_dp(self, i: int, io: int, dp: list[DP]) -> None:
		record = self.records[i]
		# observed parent
		orphan_gt = record.unphased(io)
		
		for h in range(self.NH()**2):		# hidden state
			E_all = self.emission_probability(i, h, orphan_gt)
			
			for prev_h in self.prev_h_table[h]:
				To = self.orphan_transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (To + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int], io: int) -> None:
		for i in range(self.M()):
			phased_gt = self.compute_phased_gt_by_refhaps(hs[i], i)
			self.records[i].geno[io] = phased_gt | 4
	
	def impute(self, io: int) -> None:
		L = self.NH()**2			# ハプロタイプの状態数
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = self.initialize_dp(io, L, self.M())
		for i in range(1, self.M()):
			self.update_dp(i, io, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs, io)
