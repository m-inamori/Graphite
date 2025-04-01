from __future__ import annotations

# coding: utf-8
# ParentProgenyImputer.py
# 親がphasedのときに後代をimputeする

from math import log
from typing import List, Tuple

from VCFHMM import *
from Genotype import Genotype


#################### ParentProgenyImputer ####################

MIN_PROB = -1e300

class ParentProgenyImputer(VCFHMM[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord],
							ref_haps: list[list[int]],
							is_mat_imputed: bool, map_: Map, w: float) -> None:
		VCFHMM.__init__(self, records, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_imputed = is_mat_imputed
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
	
	def num_progenies(self) -> int:
		return len(self.records[0].v) - 11
	
	# ハプロタイプを決めたときのGenotype
	def gt_by_haplotypes(self, hc1: int, hc2: int,
								mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	# i: record index, j: progeny index
	def emission_probability(self, i: int, h: int, op: int, ocs: list[int],
												phased_parent_gt: int) -> float:
		N = self.num_progenies()
		hp = h >> (N*2)
		hp2, hp1 = divmod(hp, self.NH())
		non_phased_parent_gt = (self.ref_haps[hp1][i] |
								(self.ref_haps[hp2][i] << 1))
		mat_gt = (phased_parent_gt if self.is_mat_imputed
										else non_phased_parent_gt)
		pat_gt = (non_phased_parent_gt if self.is_mat_imputed
										else phased_parent_gt)
		Ep = self.E[non_phased_parent_gt][op]	# parent emission
		# 後代の排出確率を計算
		Ec = 0.0
		for j in range(N):
			hc2, hc1 = divmod((h >> (j * 2)) & 3, 2)
			gtc = self.gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
			Ec += self.E[gtc][ocs[j]]
		return Ep + Ec
	
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
	
	def initialize_dp(self) -> list[DP]:
		M = len(self.records)	# マーカー数
		L = self.NH()**2 << (2*self.num_progenies())
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		phased_col = 9 if self.is_mat_imputed else 10
		non_phased_col = 10 if self.is_mat_imputed else 9
		record = self.records[0]
		phased_parent_gt = Genotype.phased_gt_to_int(record.v[phased_col])
		op = Genotype.gt_to_int(record.v[non_phased_col])	# observed parent
		# observed progs
		ocs = [ Genotype.gt_to_int(gt) for gt in record.v[11:] ]
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h, op, ocs, phased_parent_gt)
			dp[0][h] = (E_all, h)
		return dp
	
	# hidden stateに対して、可能な前のhidden stateを集めておく
	def collect_possible_previous_hidden_states(self) -> list[list[int]]:
		N = self.num_progenies()
		L = self.NH()**2 << (2*N)
		Lc = 1 << (2*N)				# 子のハプロタイプの状態数
		prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
		for h in range(L):		# hidden state
			hc = h & (Lc - 1)
			hp = h >> (N*2)
			hp2, hp1 = divmod(hp, self.NH())
			# 両側乗り越えることはないとする
			# non-phasedの親のあり得る前の状態
			prev_hps = []
			for h1 in range(self.NH()):
				prev_h1 = h1 + hp2 * self.NH()
				prev_hps.append(prev_h1)
			for h2 in range(self.NH()):
				if h2 == hp2:	# for duplication
					continue
				prev_h1 = hp1 + h2 * self.NH()
				prev_hps.append(prev_h1)
			
			# 後代のあり得る前の状態
			prev_hcs = []
			for hc1 in range(1 << (N*2)):
				t = hc ^ hc1		# 乗り換えたかをビットで表す
				# 乗り換えは合わせて1回まで許容できる
				if sum(1 for k in range(N*2) if ((t >> k) & 1) == 1) <= 1:
					prev_hcs.append(hc1)
			
			for prev_hp, prev_hc in product(prev_hps, prev_hcs):
				if prev_hp == hp or prev_hc == hc:
					prev_h = (prev_hp << (N*2)) | prev_hc
					prev_h_table[h].append(prev_h)
		
		return prev_h_table
	
	def compute_non_phased_parent_gt(self, h: int, i: int) -> int:
		N = self.num_progenies()
		hp = h >> (N*2)
		hp2, hp1 = divmod(hp, self.NH())
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		N = self.num_progenies()
		L = self.NH()**2 << (2*N)
		Lc = 1 << (2*N)				# 子のハプロタイプの状態数
		record = self.records[i]
		phased_col = 9 if self.is_mat_imputed else 10
		non_phased_col = 10 if self.is_mat_imputed else 9
		cc = self.Cc[i-1]	# 遷移確率
		cp = self.Cp[i-1]	# 親の遷移確率
		phased_parent_gt = Genotype.phased_gt_to_int(record.v[phased_col])
		op = Genotype.gt_to_int(record.v[non_phased_col])	# observed parent
		# observed progs
		ocs = [ Genotype.gt_to_int(gt) for gt in record.v[11:] ]
		
		for h in range(L):		# hidden state
			E_all = self.emission_probability(i, h, op, ocs, phased_parent_gt)
			
			hc = h & (Lc - 1)
			hp2, hp1 = divmod(h >> (N*2), self.NH())
			
			for prev_h in self.prev_h_table[h]:
				t1 = prev_h ^ hc	# 乗り換えしていないかしたかをbitで
				# 遷移確率 0なら(1-c)、1ならcを掛ける
				Tc = 0.0
				for j in range(N*2):
					Tc += log(cc if ((t1 >> j) & 1) == 1 else 1.0 - cc)
				# 親側の遷移
				prev_hp2, prev_hp1 = divmod(prev_h >> (N * 2), self.NH())
				Tp = (log(cp if prev_hp1 != hp1 else 1.0 - cp) +
					  log(cp if prev_hp2 != hp2 else 1.0 - cp))
				
				prob = dp[i-1][prev_h][0] + (Tc + Tp + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		M = len(self.records)
		N = self.num_progenies()
		for i in range(M):
			record = self.records[i]
			if self.is_mat_imputed:
				pat_gt = self.compute_non_phased_parent_gt(hs[i], i)
				record.set_GT(1, Genotype.int_to_phased_gt(pat_gt))
				mat_gt = Genotype.phased_gt_to_int(record.v[9])
			else:
				mat_gt = self.compute_non_phased_parent_gt(hs[i], i)
				record.set_GT(0, Genotype.int_to_phased_gt(mat_gt))
				pat_gt = Genotype.phased_gt_to_int(record.v[10])
			for j in range(N):	# 個々の後代
				hc2, hc1 = divmod((hs[i] >> (j * 2)) & 3, 2)
				gtc_int = self.gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
				record.set_GT(j+2, Genotype.int_to_phased_gt(gtc_int))
	
	def impute(self) -> None:
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		M = len(self.records)
		dp = self.initialize_dp()
		for i in range(1, M):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
