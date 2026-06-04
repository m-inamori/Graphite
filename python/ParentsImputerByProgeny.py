from __future__ import annotations

# coding: utf-8
# ParentsImputerByProgeny.py
# 後代がimputedで片親がknownな家系の補完
# 後代が複数だと難しいので一つでimputedとする

from math import log
from typing import List, Tuple

from VCFFamily import VCFFamilyRecord
from VCFHMM import *
from Genotype import Genotype


#################### ParentsImputerByProgeny ####################

MIN_PROB = -1e300

class ParentsImputerByProgeny(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
								ref_haps_mat: list[list[int]],
								ref_haps_pat: list[list[int]],
								map_: Map, w: float):
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps_mat = ref_haps_mat
		self.ref_haps_pat = ref_haps_pat
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
		return len(self.ref_haps_mat)
	
	def num_states(self) -> int:
		return self.NH()**2 << 3
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc1: int, hc2: int,
								mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	def decode_state(self, h: int) -> tuple[int, int, int, int, int]:
		hw = h & 1			# 後代のどちらに母親から渡ったか
		hc1 = (h >> 1) & 1	# 母親のどちらから後代に渡ったか
		hc2 = (h >> 2) & 1	# 父親のどちらから後代に渡ったか
		hp2, hp1 = divmod(h >> 3, self.NH())	# 親の後代に渡ってこない方
		return (hw, hc1, hc2, hp1, hp2)
	
	def compute_parent_gt(self, i: int, h: int) -> tuple[int, int]:
		def parent_gt(hc: int, hp: int, prog_a: int, b: bool) -> int:
			ref_haps = self.ref_haps_mat if b else self.ref_haps_pat
			if hc == 0:
				return prog_a | (ref_haps[hp][i] << 1)
			else:
				return ref_haps[hp][i] | (prog_a << 1)
		
		hw, hc1, hc2, hp1, hp2 = self.decode_state(h)
		prog_allele1 = self.records[i].get_allele(2, hw)	# 親から渡ったアレル
		prog_allele2 = self.records[i].get_allele(2, 1-hw)	# 親から渡ったアレル
		return (parent_gt(hc1, hp1, prog_allele1, True),
				parent_gt(hc2, hp2, prog_allele2, False))
	
	# i: record index
	def emission_probability(self, i: int, h: int) -> float:
		mat_phased_gt, pat_phased_gt = self.compute_parent_gt(i, h)
		mat_gt = self.records[i].unphased_mat()
		pat_gt = self.records[i].unphased_pat()
		return self.E[mat_phased_gt][mat_gt] + self.E[pat_phased_gt][pat_gt]
	
	def transition_probability(self, i: int, prev_h: int, h: int) -> float:
		cc = self.Cc[i-1]	# 後代の遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		hw, hc1, hc2, hp1, hp2 = self.decode_state(h)
		_, prev_hc1, prev_hc2, prev_hp1, prev_hp2 = self.decode_state(prev_h)
		
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
				log(cc if prev_hc2 != hc2 else 1.0 - cc) +
				log(cp if prev_hp1 != hp1 else 1.0 - cp) +
				log(cp if prev_hp2 != hp2 else 1.0 - cp))
	
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
			hw, hc1, hc2, hp1, hp2 = self.decode_state(h)
			# 複数が同時に乗り換えることはないとする
			for prev_h in range(hw, L, 2):
				_, prev_hc1, prev_hc2, prev_hp1, prev_hp2 = \
											self.decode_state(prev_h)
				counter = (int(prev_hc1 != hc1) +
						   int(prev_hc2 != hc2) +
						   int(prev_hp1 != hp1) +
						   int(prev_hp2 != hp2))
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
		for i in range(M):
			record = self.records[i]
			mat_gt, pat_gt = self.compute_parent_gt(i, hs[i])
			record.geno[0] = mat_gt | 4
			record.geno[1] = pat_gt | 4
	
	def impute(self) -> None:
		# DP
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
