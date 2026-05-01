from __future__ import annotations

# coding: utf-8
# ParentNoRefImputer.py
# 両親ともリファレンスにないときに片親をimputeする
# phasedの遺伝子型を考慮する

from itertools import product
from math import log
from typing import List, Tuple

from VCFGeno import *
from VCFFamily import VCFFamilyRecord
from VCFHMMRef import *
from Genotype import Genotype


#################### ParentNoRefImputer ####################

MIN_PROB = -1e300

class ParentNoRefImputer(VCFHMMRef[VCFFamilyRecord]):
	DP = List[Tuple[float, int]]	# (log of probability, prev h)
	
	def __init__(self, records: list[VCFFamilyRecord], map_: Map,
						is_mat: bool, ref_vcf: VCFGeno, w: float) -> None:
		VCFHMMRef.__init__(self, records, map_, )
		self.records: list[VCFFamilyRecord] = records
		self.is_mat: bool = is_mat		# どちらを親をimputeするか
		self.ref_haps = ref_vcf.create_ref_haps()
		self.prev_h_table = self.collect_possible_previous_hidden_states()
		self.Epc = self.calc_Epc(w)
		# crossover values
		# 親は直接の後代でないので、乗り換え確率を高くする
		K = 5						# 親の遷移確率の子に対する倍率
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
	
	# 両親と後代のphased genotypeからの後代の排出確率
	def calc_Epc(self, w: float) -> list[list[list[float]]]:
		def divide_into_alleles(gt: int) -> list[int]:
			if gt == 4:		# N/A
				return [0, 1]
			else:
				return [gt&1, gt>>1]
		
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3, N/A: 4
		# 0|1 <=> 1|0は0.1の確率で間違えるとする
		E1 = [ [ [0.0]*5 for _ in range(5) ] for _ in range(5) ]
		for mat, pat in product(range(5), repeat=2):
			orig_prob = [0.0] * 5
			as1 = divide_into_alleles(mat)
			as2 = divide_into_alleles(pat)
			for a1, a2 in product(as1, as2):
				gt = a1 | (a2 << 1)
				if gt == 1:
					orig_prob[1] += 0.25 * 0.9
					orig_prob[2] += 0.25 * 0.1
				elif gt == 2:
					orig_prob[1] += 0.25 * 0.1
					orig_prob[2] += 0.25 * 0.9
				else:
					orig_prob[gt] += 0.25
			# エラーを入れる
			for gt in range(4):
				for error_gt in range(5):
					if error_gt == 4:
						E1[mat][pat][error_gt] = w
					elif error_gt == gt:
						E1[mat][pat][error_gt] += orig_prob[gt] * (1.0 - w*2)
					else:
						E1[mat][pat][error_gt] += orig_prob[gt] * (w/3)
		
		E = [ [ [ log(p) for p in v ] for v in w ] for w in E1 ]
		return E
	
	def NH(self) -> int:
		return len(self.ref_haps)
	
	def M(self) -> int:
		return len(self.ref_haps[0])
	
	def num_progenies(self) -> int:
		return len(self.records[0].geno) - 2
	
	def compute_phased_gt_by_refhaps(self, hp: int, i: int) -> int:
		hp2, hp1 = divmod(hp, self.NH())
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def progs_emission_probability(self, i: int,
								mat_gt: int, pat_gt: int) -> float:
		record = self.records[i]
		Ec = 0.0
		for j in range(self.num_progenies()):
			gt = record.geno[j+2]
			oc = gt & 3 if gt != Genotype.NA else 4
			Ec += self.Epc[mat_gt&3][pat_gt&3][oc]
		return Ec
	
	def mat_emission_probability(self, i: int, h: int,
								mat_gt: int, pat_gt: int) -> float:
		record = self.records[i]
		phased_mat_gt = self.compute_phased_gt_by_refhaps(h, i)
		mat_gt1 = mat_gt & 3 if mat_gt != Genotype.NA else 4
		Ep = self.E[phased_mat_gt][mat_gt1]		# mat emission
		# 後代の排出確率を計算
		Ec = self.progs_emission_probability(i, phased_mat_gt, pat_gt)
		return Ep + Ec
	
	def pat_emission_probability(self, i: int, h: int,
								mat_gt: int, pat_gt: int) -> float:
		record = self.records[i]
		phased_pat_gt = self.compute_phased_gt_by_refhaps(h, i)
		pat_gt1 = pat_gt & 3 if pat_gt != Genotype.NA else 4
		Ep = self.E[phased_pat_gt][pat_gt1]		# pat emission
		# 後代の排出確率を計算
		Ec = self.progs_emission_probability(i, mat_gt, phased_pat_gt)
		return Ep + Ec
	
	# mat_gt/pat_gt: not imputed genotype
	def emission_probability(self, i: int, h: int,
								mat_gt: int, pat_gt: int) -> float:
		if self.is_mat:
			return self.mat_emission_probability(i, h, mat_gt, pat_gt)
		else:
			return self.pat_emission_probability(i, h, mat_gt, pat_gt)
	
	def parent_transition_probability(self, i: int,
											prev_hp: int, hp: int) -> float:
		cp = self.Cp[i-1]	# 親の遷移確率
		hp2, hp1 = divmod(hp, self.NH())
		prev_hp2, prev_hp1 = divmod(prev_hp, self.NH())
		# 遷移確率 0なら(1-c)、1ならcを掛ける
		return (log(cp if prev_hp1 != hp1 else 1.0 - cp) +
				log(cp if prev_hp2 != hp2 else 1.0 - cp))
	
	def initialize_dp(self, L: int, M: int) -> list[DP]:
		dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
		record = self.records[0]
		mat_gt = record.mat_gt()
		pat_gt = record.pat_gt()
		for h in range(L):		# hidden state
			E_all = self.emission_probability(0, h, mat_gt, pat_gt)
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
	
	def update_dp(self, i: int, dp: list[DP]) -> None:
		record = self.records[i]
		# observed parent
		mat_gt = record.mat_gt()
		pat_gt = record.pat_gt()
		
		for h in range(self.NH()**2):		# hidden state
			E_all = self.emission_probability(i, h, mat_gt, pat_gt)
			
			hp2, hp1 = divmod(h, self.NH())
			
			for prev_h in self.prev_h_table[h]:
				Tp = self.parent_transition_probability(i, prev_h, h)
				prob = dp[i-1][prev_h][0] + (Tp + E_all)
				dp[i][h] = max(dp[i][h], (prob, prev_h))
	
	def update_genotypes(self, hs: list[int]) -> None:
		j = 0 if self.is_mat else 1
		for i in range(len(self.ref_haps[0])):
			phased_gt = self.compute_phased_gt_by_refhaps(hs[i], i)
			self.records[i].geno[j] = phased_gt | 4
	
	def impute(self) -> None:
		L = self.NH()**2			# 親のハプロタイプの状態数
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = self.initialize_dp(L, self.M())
		for i in range(1, self.M()):
			self.update_dp(i, dp)
		
		hs = self.trace_back(dp)
		self.update_genotypes(hs)
