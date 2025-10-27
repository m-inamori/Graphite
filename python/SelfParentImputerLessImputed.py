from __future__ import annotations

# coding: utf-8
# SelfParentImputerLessImputed.py
# 自殖で後代にもimputedなサンプルがないとき親をimputeする

from math import log
from typing import List, Tuple

from GenoRecord import GenoRecord
from VCFHMM import *
from Genotype import Genotype


#################### SelfParentImputerLessImputed ####################

MIN_PROB = -1e300

class SelfParentImputerLessImputed(VCFHMM[GenoRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[GenoRecord],
							ref_haps: list[list[int]],
							map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[GenoRecord] = records
		self.ref_haps = ref_haps
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		self.Epc = self.calc_Epc(w)
		# crossover values
		K = 5		# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	# 親と後代のnon-phased genotypeからの後代の排出確率
	def calc_Epc(self, w: float) -> list[list[float]]:
		E1: list[list[float]] = [
				[ 1.0-w*2,  w/2,       w/2,    w ],
				[ 0.25-w/4, 0.5-w/2, 0.25-w/4, w ],
				[ 0.25-w/4, 0.5-w/2, 0.25-w/4, w ],
				[ w/2,      w/2,     1.0-w*2,  w ]
		]
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def M(self) -> int:
		return len(self.ref_haps[0])
	
	def num_states(self) -> int:
		return self.NH()**2
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 1
	
	def parent_genotype(self, h: int, i: int) -> int:
		h2, h1 = divmod(h, self.NH())
		return self.ref_haps[h1][i] | (self.ref_haps[h2][i] << 1)
	
	def progs_emission_probability(self, ocs: list[int],
												parent_gt: int) -> float:
		Ec = 0.0
		for oc in ocs:
			Ec += self.Epc[parent_gt][oc]
		return Ec
	
	def parent_emission_probability(self, i: int, h: int,
											parent_gt: int) -> float:
		record = self.records[i]
		phased_parent_gt = self.parent_genotype(h, i)
		return self.E[phased_parent_gt][parent_gt]
	
	# i: record index
	def emission_probability(self, i: int, h: int,
									op: int, ocs: list[int]) -> float:
		gt_parent = self.parent_genotype(h, i)
		Ep = self.E[gt_parent][op]	# parent emission
		# progenies emission
		Ec = self.progs_emission_probability(ocs, gt_parent)
		return Ep + Ec
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.num_states()
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		op = record.unphased(0)		# observed parent
		# observed progenies
		ocs = [ record.unphased(j+1) for j in range(self.num_progenies()) ]
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h, op, ocs)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		L = self.num_states()
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			hp1, hp2 = divmod(h, self.NH())
			# 両側乗り換えることはないとする
			# non-phasedの親のあり得る前の状態
			for prev_h in range(L):
				prev_hp1, prev_hp2 = divmod(prev_h, self.NH())
				if hp1 == prev_hp1 or hp2 == prev_hp2:
					prev_h_table[h].append(prev_h)
		
		return prev_h_table
	
	def transition_probability(self, i: int, h: int, prev_h: int) -> float:
		hp1, hp2 = divmod(h, self.NH())
		prev_hp1, prev_hp2 = divmod(prev_h, self.NH())
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		cp = self.Cp[i-1]	# 親の遷移確率
		return (log(cp if hp1 != prev_hp1 else 1.0 - cp) +
				log(cp if hp2 != prev_hp2 else 1.0 - cp))
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.num_states()
		record = self.records[i]
		op = record.unphased(0)		# observed parent
		# observed progs
		ocs = [ record.unphased(i+1) for i in range(self.num_progenies()) ]
		
		for h in range(L):		# hidden state
			E_all = self.emission_probability(i, h, op, ocs)
			
			for prev_h in self.prev_h_table[h]:
				T_all = self.transition_probability(i, h, prev_h)
				
				prob = dp[i-1][prev_h][0] + (T_all + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		for i in range(M):
			record = self.records[i]
			parent_gt = self.parent_genotype(hs[i], i)
			record.geno[0] = parent_gt | 4
	
	def impute(self) -> None:
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
