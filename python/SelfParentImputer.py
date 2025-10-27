from __future__ import annotations

# coding: utf-8
# SelfParentImputer.py
# 親が自殖で後代がimputedなとき親をimputeする

from math import log
from typing import List, Tuple

from GenoRecord import GenoRecord
from VCFHMM import *
from Genotype import Genotype


#################### SelfParentImputer ####################

MIN_PROB = -1e300

class SelfParentImputer(VCFHMM[GenoRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[GenoRecord],
							ref_haps: list[list[int]],
							iprog: int, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[GenoRecord] = records
		self.ref_haps = ref_haps
		self.ic = iprog		# 何番目の後代か
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		self.Epc = self.calc_Epc(w)
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(self.records, self.records[1:]) ]
		K = 5		# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	# 親と後代のnon-phased genotypeからの後代の排出確率
	def calc_Epc(self, w: float) -> list[list[float]]:
		E1: list[list[float]] = [
				[ 1.0-w*2,	w/2,	   w/2,    w ],
				[ 0.25-w/4, 0.5-w/2, 0.25-w/4, w ],
				[ 0.25-w/4, 0.5-w/2, 0.25-w/4, w ],
				[ w/2,		w/2,	 1.0-w*2,  w ]
		]
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 9
	
	def decode_state(self, h: int) -> tuple[int, int, int, int]:
		# 1bit目: 後代の左は親のどちらのハプロタイプ由来か
		# 2bit目: 後代の右は親のどちらのハプロタイプ由来か
		# 3bit以降: 親のハプロタイプの由来
		# 親のハプロタイプは少なくとも一方は後代と同じでないといけない
		hp, hc = divmod(h, 4)
		if hp < self.NH() * 2 + 4:
			h2, h1 = divmod(hp, self.NH() + 2)
		else:
			h2, h1 = divmod(hp - self.NH() * 2, 2)
		return (h1, h2, hc & 1, hc >> 1)
	
	def num_states(self) -> int:
		return (self.NH() + 1) * 16
	
	def allele(self, h: int, i: int, prog_gt: int) -> int:
		# h: 後代1, 後代2, リファレンスハプロタイプ
		if h < 2:
			return (prog_gt >> h) & 1
		else:
			return self.ref_haps[h-2][i]
	
	def parent_genotype(self, h: int, i: int, prog_gt: int) -> int:
		hp1, hp2, hc1, hc2 = self.decode_state(h)
		a1 = self.allele(hp1, i, prog_gt)
		a2 = self.allele(hp2, i, prog_gt)
		return a1 | (a2 << 1)
	
	def progeny_genotype(self, h: int, i: int, prog_gt: int) -> int:
		hp1, hp2, hc1, hc2 = self.decode_state(h)
		gt_parent = self.parent_genotype(h, i, prog_gt)
		return ((gt_parent >> hc1) & 1) | ((gt_parent >> hc2) & 1)
	
	# phasedなGenotypeに対してphasedなGenotypeが排出される確率
	def phased_emission_probability(self, h: int, i: int,
												prog_gt: int) -> float:
		phased_prog_gt = self.progeny_genotype(h, i, prog_gt)
		d = phased_prog_gt ^ prog_gt
		return log((0.99 if (d & 1) == 0 else 0.01) *
					(0.99 if (d >> 1) == 0 else 0.01))
	
	def non_phased_emission_probability(self, ocs: list[int],
												parent_gt: int) -> float:
		return sum(self.Epc[parent_gt][oc] for oc in ocs)
	
	# i: record index
	def emission_probability(self, h: int, i: int, op: int,
									prog_gt: int, ocs: list[int]) -> float:
		# imputedな後代も排出確率を計算する
		parent_gt = self.parent_genotype(h, i, prog_gt)
		Ep = self.E[parent_gt][op]
		Ec = self.phased_emission_probability(h, i, prog_gt)
		Ecs = self.non_phased_emission_probability(ocs, parent_gt)
		return Ep + Ec + Ecs
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.num_states()
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		op = record.unphased(0)		# observed parent
		prog_gt = record.geno[self.ic+1] & 3	# phased progeny
		ocs = [ record.unphased(j+1) for j in range(self.num_progenies())
														if j != self.ic ]
		for h in range(L):		# hidden state
			E_all = self.emission_probability(h, 0, op, prog_gt, ocs)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		# 親の両側ともリファレンスハプロタイプということはない
		L = self.num_states()
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			hp1, hp2, hc1, hc2 = self.decode_state(h)
			# 両側乗り換えることはないとする
			# non-phasedの親のあり得る前の状態
			for prev_h in range(L):
				prev_hp1, prev_hp2, prev_hc1, prev_hc2 = \
												self.decode_state(prev_h)
				diff_counter = ((1 if hp1 != prev_hp1 else 0) +
								(1 if hp2 != prev_hp2 else 0) +
								(1 if hc1 != prev_hc1 else 0) +
								(1 if hc2 != prev_hc2 else 0))
				if diff_counter <= 1:
					prev_h_table[h].append(prev_h)
		
		return prev_h_table
	
	def transition_probability(self, h: int, prev_h: int,
										cc: float, cp: float) -> float:
		hp1, hp2, hc1, hc2 = self.decode_state(h)
		prev_hp1, prev_hp2, prev_hc1, prev_hc2 = self.decode_state(prev_h)
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cp if hp1 != prev_hp1 else 1.0 - cp) +
				log(cp if hp2 != prev_hp2 else 1.0 - cp) +
				log(cc if hc1 != prev_hc1 else 1.0 - cc) +
				log(cc if hc2 != prev_hc2 else 1.0 - cc))
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.num_states()
		record = self.records[i]
		cc = self.Cc[i-1]	# 遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		op = record.unphased(0)		# observed parent
		prog_gt = record.geno[self.ic+1] & 3	# phased progeny
		ocs = [ record.unphased(j+1) for j in range(self.num_progenies())
															if j != self.ic ]
		
		for h in range(L):		# hidden state
			E_all = self.emission_probability(h, i, op, prog_gt, ocs)
			
			for prev_h in self.prev_h_table[h]:
				T_all = self.transition_probability(h, prev_h, cc, cp)
				
				prob = dp[i-1][prev_h][0] + (T_all + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		for i in range(M):
			record = self.records[i]
			prog_gt = record.geno[self.ic+1] & 3
			parent_gt = self.parent_genotype(hs[i], i, prog_gt)
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
