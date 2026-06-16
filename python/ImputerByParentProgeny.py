from __future__ import annotations

# coding: utf-8
# ImputerByParentProgeny.py
# 後代がimputedで片親がimputedもう片親がknownな家系の補完
# 残りの後代も補完する

from math import log
from typing import List, Tuple

from VCFFamily import VCFFamilyRecord
from VCFHMM import *
from Genotype import Genotype


#################### ImputerByParentProgeny ####################

MIN_PROB = -1e300

class ImputerByParentProgeny(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
								ref_haps: list[list[int]],
								should_impute_mat: bool, map_: Map, w: float):
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.should_impute_mat = should_impute_mat
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
		K = 5						# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def should_imputed_index(self) -> int:
		return int(not self.should_impute_mat)
	
	def num_should_impute_progenies(self) -> int:
		return len(self.records[0].geno) - 3
	
	def num_states(self) -> int:
		return len(self.ref_haps) << (self.num_should_impute_progenies()*2+3)
	
	def decode(self, h: int) -> tuple[int, int, int, int, list[int]]:
		hw = h & 1
		hp1 = (h >> 1) & 1
		hp21 = (h >> 2) & 1
		np = self.num_should_impute_progenies()
		hc = [ (h >> i) & 1 for i in range(3, np*2+3) ]
		hp22 = h >> (np*2 + 3)
		return (hw, hp1, hp21, hp22, hc)
	
	# ハプロタイプを決めたときの親のGenotype
	def parent_gt_by_haplotypes(self, i: int, hw: int, hp1: int,
											hp2: int, prog_gt: int) -> int:
		if hp1 == 0:
			return ((prog_gt >> (1 - hw)) & 1) | (self.ref_haps[hp2][i] << 1)
		else:
			return (self.ref_haps[hp2][i] << 1) | ((prog_gt >> (1 - hw)) & 1)
	
	# ハプロタイプを決めたときの後代のGenotype
	def progeny_gt_by_haplotypes(self, j: int, hc: list[int],
										mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc[j*2]) & 1) | ((pat_gt >> hc[j*2+1]) & 1)
	
	def compute_gts(self, i: int, h: int) -> list[int]:
		gts = [0] * (self.num_should_impute_progenies() + 1)
		hw, hp1, hp21, hp22, hc = self.decode(h)
		record = self.records[i]
		prog_gt = record.geno[2]
		parent_gt = self.parent_gt_by_haplotypes(i, hw, hp21, hp22, prog_gt)
		gts[0] = parent_gt
		mat_gt = parent_gt if self.should_impute_mat else record.geno[0]
		pat_gt = record.geno[1] if self.should_impute_mat else parent_gt
		for j in range(self.num_should_impute_progenies()):
			gts[j+1] = self.progeny_gt_by_haplotypes(j, hc, mat_gt, pat_gt)
		return gts
	
	# i: record index
	def emission_probability(self, i: int, h: int) -> float:
		record = self.records[i]
		phased_gts = self.compute_gts(i, h)
		hw, hp1, hp21, hp22, hc = self.decode(h)
		index = self.should_imputed_index()
		# imputedな親からimputedな後代へ渡っているアレルは一致しているか
		parent_a = (record.geno[1-index] >> hp1) & 1
		progeny_a = (record.geno[2] >> hw) & 1
		E = log(0.99) if parent_a == progeny_a else log(0.01)
		parent_gt = record.unphased(index)
		E += self.E[phased_gts[0]][parent_gt]
		for j in range(self.num_should_impute_progenies()):
			gt = record.unphased(j+3)
			E += self.E[phased_gts[j+1]][gt]
		return E
	
	def transition_probability(self, i: int, prev_h: int, h: int) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		_, prev_hp1, prev_hp21, prev_hp22, prev_hc = self.decode(prev_h)
		__, hp1, hp21, hp22, hc = self.decode(h)
		
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		T = (log(cc if prev_hp1 != hp1 else 1.0 - cc) +
			 log(cc if prev_hp21 != hp21 else 1.0 - cc) +
			 log(cp if prev_hp22 != hp22 else 1.0 - cp))
		for j in range(len(hc)):
			T += log(cc if prev_hc[j] != hc[j] else 1.0 - cc)
		return T
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.num_states()
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		L = self.num_states()
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			hw, hp1, hp21, hp22, hc = self.decode(h)
			# 複数が同時に乗り換えることはないとする
			for prev_h in range(hw, L, 2):
				_, prev_hp1, prev_hp21, prev_hp22, prev_hc = self.decode(prev_h)
				counter = (int(prev_hp1 != hp1) +
						   int(prev_hp21 != hp21) +
						   int(prev_hp22 != hp22))
				for j in range(len(hc)):
					counter += int(prev_hc[j] != hc[j])
				if counter <= 1:
					prev_h_table[h].append(prev_h)
		
		return prev_h_table
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.num_states()
		for h in range(L):		# hidden state
			E_all = self.emission_probability(i, h)
			for prev_h in self.prev_h_table[h]:
				Tp = self.transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (Tp + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		index = self.should_imputed_index()
		for i in range(M):
			record = self.records[i]
			gts = self.compute_gts(i, hs[i])
			record.geno[index] = gts[0] | 4
			for j in range(3, self.num_should_impute_progenies()+3):
				record.geno[j] = gts[j-2] | 4
	
	def impute(self) -> None:
		# DP
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
