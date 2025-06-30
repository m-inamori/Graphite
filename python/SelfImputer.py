from __future__ import annotations

# coding: utf-8
# SelfParentImputer.py
# 自殖で親と後代を一度にimputeする

from math import log
from typing import List, Tuple

from VCFHMM import *
from Genotype import Genotype


#################### SelfParentImputerLessImputed ####################

MIN_PROB = -1e300

class SelfParentImputerLessImputed(VCFHMM[VCFRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFRecord],
							ref_haps: list[list[int]],
							map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFRecord] = records
		self.ref_haps = ref_haps
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
		K = 5		# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def M(self) -> int:
		return len(self.ref_haps[0])
	
	def num_states(self) -> int:
		return self.NH()**2 << (self.num_progenies() * 2)
	
	def num_progenies(self) -> int:
		return len(self.records[0].v) - 10
	
	def compute_parent_phased_gt(self, hp: int, i: int) -> int:
		hp2, hp1 = divmod(hp, self.NH())
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def compute_progeny_phased_gts(self, hc: int, parent_gt: int) -> list[int]:
		return [ ((parent_gt >> (hc >> (j*2))) & 1) |
				 ((parent_gt >> (hc >> (j*2+1))) & 1)
							for j in range(self.num_progenies()) ]
	
	def progs_emission_probability(self, hc: int, ocs: list[int],
												parent_gt: int) -> float:
		phased_gts = self.compute_progeny_phased_gts(hc, parent_gt)
		Ec = 0.0
		for j in range(self.num_progenies()):
			Ec += self.E[phased_gts[j]][ocs[j]]
		return Ec
	
	def parent_emission_probability(self, i: int, h: int,
											parent_gt: int) -> float:
		record = self.records[i]
		hc, hp = divmod(h, self.NH()**2)
		phased_parent_gt = self.compute_parent_phased_gt(h, i)
		return self.E[phased_parent_gt][parent_gt]
	
	def parent_genotype(self, hp: int, i: int) -> int:
		h2, h1 = divmod(hp, self.NH())
		return self.ref_haps[h1][i] | (self.ref_haps[h2][i] << 1)
	
	# i: record index
	def emission_probability(self, i: int, h: int,
									op: int, ocs: list[int]) -> float:
		hc, hp = divmod(h, self.NH()**2)
		gt_parent = self.parent_genotype(hp, i)
		Ep = self.E[gt_parent][op]	# parent emission
		# progenies emission
		Ec = self.progs_emission_probability(hc, ocs, gt_parent)
		return Ep + Ec
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.num_states()
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		op = Genotype.gt_to_int(record.v[9])	# observed parent
		# observed progenies
		ocs = [ Genotype.gt_to_int(record.v[j+10])
							for j in range(self.num_progenies()) ]
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
	
	def parent_transition_probability(self, i: int, hp: int,
													prev_hp: int) -> float:
		hp1, hp2 = divmod(hp, self.NH())
		prev_hp1, prev_hp2 = divmod(prev_hp, self.NH())
		cp = self.Cp[i-1]	# 親の遷移確率
		return (log(cp if hp1 != prev_hp1 else 1.0 - cp) +
				log(cp if hp2 != prev_hp2 else 1.0 - cp))
	
	def progeny_transition_probability(self, i: int, hc: int,
													prev_hc: int) -> float:
		T = 0.0
		cc = self.Cc[i-1]	# 後代の遷移確率
		for j in range(self.num_progenies()*2):
			if ((hc >> j) & 1) != ((prev_hc >> j) & 1):
				T += log(cc)
			else:
				T += log(1.0 - cc)
		return T
	
	def transition_probability(self, i: int, h: int, prev_h: int) -> float:
		hc, hp = divmod(h, self.NH()**2)
		prev_hc, prev_hp = divmod(prev_h, self.NH()**2)
		Tp = self.parent_transition_probability(i, hp, prev_hp)
		Tc = self.progeny_transition_probability(i, hc, prev_hc)
		return Tp + Tc
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.num_states()
		record = self.records[i]
		op = Genotype.gt_to_int(record.v[9])	# observed parent
		# observed progs
		ocs = [ Genotype.phased_gt_to_int(record.v[i+10])
						for i in range(self.num_progenies()) ]
		
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
			hc, hp = divmod(hs[i], self.NH()**2)
			parent_gt = self.parent_genotype(hp, i)
			record.set_GT(0, Genotype.int_to_phased_gt(parent_gt))
	
	def impute(self) -> None:
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
