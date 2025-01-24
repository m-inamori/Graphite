from __future__ import annotations

# coding: utf-8
# VCFNoParentImputed.py
# 両親ともimputeされていない
# 母親を補完する

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple

from VCFFamily import *
from Map import *
from Genotype import Genotype


#################### VCFNoParentImputed ####################

DP = List[Tuple[float, int]]	# (log of probability, prev h)
MIN_PROB = -1e300

class VCFNoParentImputed(VCFBase, VCFSmallBase,
								VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
										ref_haps: list[list[int]], map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.E = self.calc_E()
		self.Epc = self.calc_Epc()
		# crossover values
		# 親は直接の後代でないので、乗り換え確率を高くする
		K = 5						# 親の遷移確率の子に対する倍率
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]

	def calc_E(self):
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E1 = [ [ 1.0-w*2, w/2,     w/2,     w ],
			   [ w/2,     1.0-w*2, w/2,     w ],
			   [ w/2,     1.0-w*2, w/2,     w ],
			   [ w/2,     w/2,     1.0-w*2, w ] ]
		
		E = [ [ log(p) for p in v ] for v in E1 ]
		return E
	
	# 両親と後代のnon-phased genotypeからの後代の排出確率
	def calc_Epc(self):
		w = 0.01		# probability of wrong
		E1 = [ [ [ 1.0-w*2,   w/2,       w/2,       w ],
				 [ 0.5-w*3/4, 0.5-w*3/4, w/2,       w ],
				 [ w/2,       1.0-w*2,   w/2,       w ],
				 [ 0.5-w*3/4, 0.5-w*3/4, w/2,       w ] ],
			   [ [ 0.5-w*3/4, 0.5-w*3/4, w/2,       w ],
			     [ 1/4-w/8,   0.5-w*3/4, 1/4-w/8,   w ],
			     [ w/2,       0.5-w*3/4, 0.5-w*3/4, w ],
			     [ 1/4+w/4,   0.5-w*3/4, 1/4+w/4,   w ] ],
			   [ [ w/2,       1.0-w*2,   w/2,       w ],
			     [ w/2,       0.5-w*3/4, 0.5-w*3/4, w ],
			     [ w/2,       w/2,       1.0-w*2,   w ],
			     [ w/2,       0.5-w*3/4, 0.5-w*3/4, w ] ] ]
		E = [ [ [ log(p) for p in v ] for v in w ] for w in E1 ]
		return E
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def dist(self, r1: VCFRecord, r2: VCFRecord) -> float:
		dist = (self.cM(r2.pos()) - self.cM(r1.pos())) / 100
		if dist > 0.0:
			return dist
		else:
			return (r2.pos() - r1.pos()) * 1e-6
	
	def impute(self):
		# mat_gt/pat_gt: not imputed genotype
		def emission_probability(i: int, h: int,
									mat_gt: int, pat_gt: int) -> float:
			hp2, hp1 = divmod(h, NH)
			record = self.records[i]
			phased_mat_gt = self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
			Ep = self.E[phased_mat_gt][mat_gt]		# mat emission
			# 後代の排出確率を計算
			Ec = 0.0
			non_phased_mat_gt = (phased_mat_gt & 1) + (phased_mat_gt >> 1)
			for j in range(N):
				oc = Genotype.gt_to_int(record.v[j+11])
				Ec += self.Epc[non_phased_mat_gt][pat_gt][oc]
			return Ep + Ec
		
		def initialize_dp(L: int, M: int) -> list[DP]:
			dp = [ [ (MIN_PROB, 0) ] * L for _ in range(M) ]
			record = self.records[0]
			mat_gt = Genotype.gt_to_int(record.v[9])
			pat_gt = Genotype.gt_to_int(record.v[10])
			for h in range(L):		# hidden state
				E_all = emission_probability(0, h, mat_gt, pat_gt)
				dp[0][h] = (E_all, h)
			return dp
		
		# hidden stateに対して、可能な前のhidden stateを集めておく
		def collect_possible_previous_hidden_states(L: int) -> list[list[int]]:
			prev_h_table: list[list[int]] = [ [] for _ in range(L) ]
			for h in range(L):		# hidden state
				hp2, hp1 = divmod(h, NH)
				# 両側乗り換えることはないとする
				# non-phasedの親のあり得る前の状態
				for h1 in range(NH):
					prev_h1 = h1 + hp2 * NH
					prev_h_table[h].append(prev_h1)
				for h2 in range(NH):
					if h2 == hp2:	# for duplication
						continue
					prev_h1 = hp1 + h2 * NH
					prev_h_table[h].append(prev_h1)
			
			return prev_h_table
		
		def update_dp(i: int, dp: list[DP]):
			record = self.records[i]
			cp = self.Cp[i-1]	# 親の遷移確率
			# observed parent
			mat_gt = Genotype.gt_to_int(record.v[9])
			pat_gt = Genotype.gt_to_int(record.v[10])
			
			for h in range(L):		# hidden state
				E_all = emission_probability(i, h, mat_gt, pat_gt)
				
				hp2, hp1 = divmod(h, NH)
				
				for prev_h in prev_h_table[h]:
					prev_hp2, prev_hp1 = divmod(prev_h, NH)
					# 遷移確率 0なら(1-c)、1ならcを掛ける
					Tp = (log(cp if prev_hp1 != hp1 else 1.0 - cp) +
						  log(cp if prev_hp2 != hp2 else 1.0 - cp))
					
					prob = dp[i-1][prev_h][0] + (Tp + E_all)
					dp[i][h] = max(dp[i][h], (prob, prev_h))
		
		def trace_back(dps: list[list[float]]) -> list[int]:
			hs = [ 0 ] * M
			h, (_, prev_h) = max(enumerate(dp[-1]), key=lambda v: v[1])
			hs[-1] = h
			hs[-2] = prev_h
			for i in range(M-2, 0, -1):
				_, prev_h = dp[i][prev_h]
				hs[i-1] = prev_h
			return hs
		
		def compute_non_phased_parent_gt(h: int, i: int) -> int:
			hp2, hp1 = divmod(h, NH)
			return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
		
		def update_genotypes(hs: list[int]):
			for i in range(M):
				phased_mat_gt = compute_non_phased_parent_gt(hs[i], i)
				self.records[i].v[9] = Genotype.int_to_phased_gt(phased_mat_gt)
		
		N = self.num_progenies()
		NH = len(self.ref_haps)		# リファレンスのハプロタイプの数
		L = NH**2					# 親のハプロタイプの状態数
		M = len(self)				# マーカー数
		
		prev_h_table = collect_possible_previous_hidden_states(L)
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = initialize_dp(L, M)
		for i in range(1, M):
			update_dp(i, dp)
		
		hs = trace_back(dp)
		update_genotypes(hs)
