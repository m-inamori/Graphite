from __future__ import annotations

# coding: utf-8
# VCFOneParentImputed.py
# 片側だけimputeされている

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from Map import *


#################### VCFOneParentImputed ####################

class VCFOneParentImputed(VCFBase, VCFSmallBase, VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_imputed = is_mat_imputed
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def impute(self):
		def phased_gt_to_int(gt: str) -> int:
			gt1 = 0 if gt[0] == '0' else 1
			gt2 = 0 if gt[2] == '0' else 1
			return gt1 | (gt2 << 1)
		
		def gt_to_int(gt: str) -> int:
			if '.' in gt:
				return 3
			gt1 = 0 if gt[0] == '0' else 1
			gt2 = 0 if gt[2] == '0' else 1
			return gt1 + gt2
		
		def int_to_phased_gt(gt_int: int) -> str:
			if gt_int == 0:
				return '0|0'
			elif gt_int == 1:
				return '1|0'
			elif gt_int == 2:
				return '0|1'
			else:
				return '1|1'
		
		# ハプロタイプを決めたときのGenotype
		def gt_by_haplotypes(hc1: int, hc2: int,
								mat_gt: int, pat_gt: int) -> int:
			return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
		
		def emission_probability(i: int, h: int, op: int, ocs: list[int],
												phased_parent_gt: int) -> float:
			hp = h >> (N*2)
			hp2, hp1 = divmod(hp, NH)
			non_phased_parent_gt = (self.ref_haps[hp1][i] |
									(self.ref_haps[hp2][i] << 1))
			mat_gt = (phased_parent_gt if self.is_mat_imputed
											else non_phased_parent_gt)
			pat_gt = (non_phased_parent_gt if self.is_mat_imputed
											else phased_parent_gt)
			Ep = E[non_phased_parent_gt][op]	# parent emission
			# 後代の排出確率を計算
			Ec = 0.0
			for j in range(N):
				hc2, hc1 = divmod((h >> (j * 2)) & 3, 2)
				gtc = gt_by_haplotypes(hc1, hc2, non_phased_parent_gt,
														phased_parent_gt)
				Ec += E[gtc][ocs[j]]
			return Ep + Ec
		
		DP = List[Tuple[float, int]]	# (log of probability, prev h)
		
		def initialize_dp(L: int, M: int) -> list[DP]:
			dp = [ [ (-1e300, 0) ] * L for _ in range(M) ]
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			record = self.records[0]
			phased_parent_gt = phased_gt_to_int(record.v[phased_col])
			op = gt_to_int(record.v[non_phased_col])		# observed parent
			# observed progs
			ocs = [ gt_to_int(gt) for gt in record.v[11:] ]
			for h in range(L):		# hidden state
				E_all = emission_probability(0, h, op, ocs, phased_parent_gt)
				dp[0][h] = (E_all, h)
			return dp
		
		def update_dp(i: int, dp: list[DP]):
			record = self.records[i]
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			cc = Cc[i-1]	# 遷移確率
			cp = Cp[i-1]	# 親の遷移確率
			phased_parent_gt = phased_gt_to_int(record.v[phased_col])
			op = gt_to_int(record.v[non_phased_col])		# observed parent
			# observed progs
			ocs = [ gt_to_int(gt) for gt in record.v[11:] ]
			
			for h in range(L):		# hidden state
				E_all = emission_probability(i, h, op, ocs, phased_parent_gt)
				
				hc = h & (Lc - 1)
				hp2, hp1 = divmod(h >> (N*2), NH)
				# 両側乗り越えることはないとする
				# non-phasedの親のあり得る前の状態
				prev_hps = []
				for h1 in range(NH):
					prev_h1 = ((h1 + hp2 * NH) << (N*2))
					prev_hps.append(prev_h1)
				for h2 in range(NH):
					if h2 == hp2:	# for duplication
						continue
					prev_h1 = ((hp1 + h2 * NH) << (N*2))
					prev_hps.append(prev_h1)
				
				# 後代のあり得る前の状態
				prev_hcs = []
				for hc1 in range(1 << (N*2)):
					if all((((hc1 >> (k*2)) & 3) ^ hc) != 3 for k in range(N)):
						prev_hcs.append(hc1)
				
				for prev_hp, prev_hc in product(prev_hps, prev_hcs):
					prev_h = prev_hp + prev_hc
					t1 = prev_h ^ hc	# 乗り換えしていないかしたかをbitで
					# 遷移確率 0なら(1-c)、1ならcを掛ける
					Tc = 0.0
					for j in range(N*2):
						Tc += log(cc if ((t1 >> j) & 1) == 1 else 1.0 - cc)
					# 親側の遷移
					prev_hp2, prev_hp1 = divmod(prev_h >> (N * 2), NH)
					Tc += log(cp if prev_hp1 != hp1 else 1.0 - cp)
					Tc += log(cp if prev_hp2 != hp2 else 1.0 - cp)
					
					dp[i][h] = max(dp[i][h],
									(dp[i-1][prev_h][0] + Tc + E_all, prev_h))
		
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
			hp = hs[i] >> (N*2)
			hp2, hp1 = divmod(hp, NH)
			return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
		
		def update_genotypes(hs: list[int]):
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			for i in range(M):
				record = self.records[i]
				non_phased_parent_gt = compute_non_phased_parent_gt(hs[i], i)
				record.v[non_phased_col] = int_to_phased_gt(
													non_phased_parent_gt)
				phased_parent_gt = phased_gt_to_int(record.v[phased_col])
				mat_gt = (phased_parent_gt if self.is_mat_imputed
												else non_phased_parent_gt)
				pat_gt = (non_phased_parent_gt if self.is_mat_imputed
												else phased_parent_gt)
				for j in range(N):	# 個々の後代
					hc2, hc1 = divmod((hs[i] >> (j * 2)) & 3, 2)
					gtc_int = gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
					record.v[j+11] = int_to_phased_gt(gtc_int)
		
		N = self.num_progenies()
		M = len(self)				# マーカー数
		K = 5						# 親の遷移確率の子に対する倍率
		NH = len(self.ref_haps)		# リファレンスのハプロタイプの数
#		M = 10	# ****** 仮
#		NH = 3	# ****** 仮
		Lp = NH**2					# 親のハプロタイプの状態数
		Lc = 1 << (2*N)				# 子のハプロタイプの状態数
		L = Lp * Lc
		
		# crossover values
		# 後代
		Cc = [ Map.Kosambi((self.cM(r2.pos()) - self.cM(r1.pos())) / 100)
							for r1, r2 in zip(self.records, self.records[1:]) ]
		# 親は直接の後代でないので、乗り換え確率を高くする
		Cp = [ Map.Kosambi((self.cM(r2.pos()) - self.cM(r1.pos())) / 100 * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
		
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E1 = [ [ 1.0-3*w, w,       w,       w ],
			   [ w,       1.0-3*w, w,       w ],
			   [ w,       1.0-3*w, w,       w ],
			   [ w,       w,       1.0-3*w, w ] ]
		
		E = [ [ log(p) for p in v ] for v in E1 ]
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = initialize_dp(L, M)
		for i in range(1, M):
			update_dp(i, dp)
		
		hs = trace_back(dp)
		update_genotypes(hs)
