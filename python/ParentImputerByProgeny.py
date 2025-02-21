from __future__ import annotations

# coding: utf-8
# ParentImputerByProgeny.py
# 後代がimputedで片親がknownな家系の補完
# 後代が複数だと難しいので一つでimputedとする

from math import log
from typing import List, Tuple

from VCFHMM import *
from Genotype import Genotype


#################### ParentImputerByProgeny ####################

MIN_PROB = -1e300

class ParentImputerByProgeny(VCFHMM[VCFRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFRecord], ref_haps: list[list[int]],
							is_mat_known: bool, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_known = is_mat_known
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
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc1: int, hc2: int,
								mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	def compute_parent_gt(self, i: int, h: int) -> int:
		hw = h & 1			# 後代のどちらに親から渡ったか
		hc = (h >> 1) & 1	# 親のどちらから後代に渡ったか
		hp = h >> 2			# 後代に渡っていない方
		prog_allele = int(self.records[i].v[10][hw*2])	# 親から渡ったアレル
		if hc == 0:
			return prog_allele | (self.ref_haps[hp][i] << 1)
		else:
			return self.ref_haps[hp][i] | (prog_allele << 1)
	
	# i: record index
	def emission_probability(self, i: int, h: int) -> float:
		parent_phased_gt = self.compute_parent_gt(i, h)
		parent_gt = Genotype.gt_to_int(self.records[i].v[9])
		return self.E[parent_phased_gt][parent_gt]	# parent emission
	
	def transition_probability(self, i: int, prev_h: int, h: int) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		hc = (h >> 1) & 1
		hp = h >> 2
		prev_hc = (prev_h >> 1) & 1
		prev_hp = prev_h >> 2
		
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc != hc else 1.0 - cc) +
				log(cp if prev_hp != hp else 1.0 - cp))
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.NH() << 2
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		L = self.NH() << 2
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			# 後代のどちらに親から渡ってくるか
			# これは不変
			hw = h & 1
			hc = (h >> 1) & 1	# 親のどちらが渡ってくるか
			hp = h >> 2			# 親の後代に渡ってこない方
			# 複数が同時に乗り換えることはないとする
			for prev_h in range(hw, L, 2):
				prev_hc = (prev_h >> 1) & 1
				prev_hp = prev_h >> 2
				if prev_hp == hp or prev_hc == hc:
					prev_h_table[h].append(prev_h)
		
		return prev_h_table
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.NH() << 2
		for h in range(L):		# hidden state
			E_all = self.emission_probability(i, h)
			for prev_h in self.prev_h_table[h]:
				Tp = self.transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (Tp + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		# 後代のどちらに渡ったかでGenotypeが入れ替わる可能性がある
		is_swapped = self.is_mat_known ^ (hs[0] & 1) == 0
		for i in range(M):
			record = self.records[i]
			parent_gt = self.compute_parent_gt(i, hs[i])
			record.v[9] = Genotype.int_to_phased_gt(parent_gt)
			if is_swapped:
				gt = record.v[10]
				record.v[10] = gt[2] + '|' + gt[0] + gt[3:]
	
	def impute(self) -> None:
		# DP
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
