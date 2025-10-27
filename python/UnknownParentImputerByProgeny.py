from __future__ import annotations

# coding: utf-8
# ProgenyImputerByProgeny.py
# 後代が一つimputedで片親がunknownな家系の補完
# unknownな親とnon imputedな後代を一つ補完する

from math import log
from typing import List, Tuple

from VCFFamily import VCFFamilyRecord
from VCFHMM import *
from Genotype import Genotype


#################### UnknownParentImputerByProgeny ####################

MIN_PROB = -1e300

class ParentImputerByProgeny(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
							ref_haps: list[list[int]],
							is_mat_known: bool, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_known = is_mat_known
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		# crossover values
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(records, records[1:]) ]
		# 後代から親に行ってまた後代へ
		self.Cc2 = [ Map.Kosambi(self.dist(r1, r2) * 2)
							for r1, r2 in zip(records, records[1:]) ]
		K = 5						# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def num_progenies(self) -> int:
		return 1
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc1: int, hc2: int,
								mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	def compute_progeny_gt(self, i: int, h: int) -> int:
		record = self.records[i]
		hc1 = h & 1			# imputedの親のどちらが後代に渡ったか
		hc2 = (h >> 1) & 1	# imputedな後代とref_hapsのどちらが後代に渡ったか
		hp = h >> 2			# ref_hapsからどれが来たか
		a1 = record.get_allele(0, hc1)	# imputedな親から渡ったアレル
		k = 2 if self.is_mat_known else 0	# unknownな親から渡った側
		a2 = record.get_allele(1, k) if hc2 == 0 else self.ref_haps[i][hp]
		return a1 | (a2 << 1) if self.is_mat_known else a2 | (a1 << 1)
	
	# i: record index
	def emission_probability(self, i: int, h: int) -> float:
		hc1 = h & 1			# imputedの親のどちらが後代に渡ったか
		hc2 = (h >> 1) & 1	# imputedな後代とref_hapsのどちらが後代に渡ったか
		hp = h >> 2			# ref_hapsからどれが来たか
		parent_phased_gt = self.compute_parent_gt(i, h)
		parent_gt = self.records[i].unphased(self.non_phased_col())
		return self.E[parent_phased_gt][parent_gt]	# parent emission
	
	def transition_probability(self, i: int, prev_h: int, h: int) -> float:
		cc = self.Cc[i-1]		# 後代の遷移確率
		cc2 = self.Cc2[i-1]		# 後代⇒親⇒後代の遷移確率
		cp = self.Cp[i-1]		# 親の遷移確率
		hc1 = h & 1
		hc2 = (h >> 1) & 1
		hp = h >> 2
		prev_hc1 = prev_h & 1
		prev_hc2 = (prev_h >> 1) & 1
		prev_hp = prev_h >> 2
		
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cc if prev_hc1 != hc1 else 1.0 - cc) +
				log(cc2 if prev_hc2 != hc2 else 1.0 - cc2) +
				log(cp if prev_hp != hp else 1.0 - cp))
	
	def phased_col(self) -> int:
		return 0 if self.is_mat_known else 1
	
	def non_phased_col(self) -> int:
		return 1 if self.is_mat_known else 0
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.NH() << 2
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		phased_parent_gt = record.geno[2] & 3
		op = record.unphased(self.non_phased_col())		# observed parent
		# observed progs
		ocs = [ record.unphased(k) for k in range(2, len(record.geno)) ]
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
			prev_hps = []
			for prev_h in range(hw, L, 2):
				prev_hc = (prev_h >> 1) & 1
				prev_hp = prev_h >> 2
				if prev_hp == hp or prev_hc == hc:
					prev_hps.append(prev_h)
		
		return prev_h_table
	
	def compute_parent_gt(self, h: int, i: int) -> int:
		N = self.num_progenies()
		hp = h >> (N*2)
		hp2, hp1 = divmod(hp, self.NH())
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		L = self.NH() << 2
		record = self.records[i]
		cc = self.Cc[i-1]	# 遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		phased_parent_gt = record.geno[self.phased_col()] & 3
		op = record.unphased(self.non_phased_col())		# observed parent
		
		for h in range(L):		# hidden state
			E_all = self.emission_probability(i, h)
			for prev_h in self.prev_h_table[h]:
				T_all = self.transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (T_all + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		# 後代のどちらに渡ったかでGenotypeが入れ替わる可能性がある
		is_swapped = self.is_mat_known ^ (hs[0] & 1) == 0
		for i in range(M):
			record = self.records[i]
			parent_gt = self.compute_parent_gt(i, hs[i])
			record.geno[0] = parent_gt | 4
			if is_swapped:
				gt = record.geno[1]
				record.geno[1] = Genotype.inverse(gt)
	
	def impute(self) -> None:
		# DP
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
