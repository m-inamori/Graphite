from __future__ import annotations

# coding: utf-8
# VCFOneParentImputed.py
# 片側だけimputeされている

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from Map import *


#################### VCFOneParentImputed ####################

class VCFOneParentImputed(VCFBase, VCFSmallBase, VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]],
						records: list[VCFFamilyRecord],
						is_mat_imputed: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
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
		# hc: hidden state of progeny
		# hp: hidden state of parent
		def gt_by_haplotype(hc: int, hp: int, phased_parent_gt) -> int:
			hap1 = hc & 1
			hap2 = hc >> 1
			if self.is_mat_imputed:
				mat = (phased_parent_gt >> hap1) & 1
				pat = (hp >> hap2) & 1
			else:
				mat = (hp >> hap1) & 1
				pat = (phased_parent_gt >> hap2) & 1
			return mat | (pat << 1)
		
		def initialize_dp(L: int, M: int) -> list[list[float]]:
			dp = [ [ 0.0 ] * L for _ in range(M + 1) ]
			for j in range(L):
				dp[0][j] = 1.0
			return dp
		
		def emission_probability(h: int, op: int, ocs: list[int],
											phased_parent_gt: int) -> float:
			hp = h >> (N*2)
			hc = h & ((1 << (N*2)) - 1)
			Ep = E[hp][op]			# parent emission
			# 後代の排出確率を計算
			Ec = 1.0
			for j in range(N):
				hc = (h >> (j * 2)) & 3
				gtc = gt_by_haplotype(hc, hp, phased_parent_gt)
				Ec *= E[gtc][ocs[j]]
			return Ep * Ec
		
		def update_dp(i: int, dp: list[list[float]]):
			record = self.records[i]
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			c = C[i]	# 遷移確率
			phased_parent_gt = phased_gt_to_int(record.v[phased_col])
			op = gt_to_int(record.v[non_phased_col])		# observed parent
			ocs = [ gt_to_int(gt) for gt in record.v[11:] ]		# observed progs
			
			for h in range(L):		# hidden state
				E_all = emission_probability(h, op, ocs, phased_parent_gt)
				
				hc = h & ((1 << (N*2)) - 1)
				for prev_h in range(L):		# previous hidden state
					t = prev_h ^ hc		# 乗り換えしていないかしたかをbitで
					# 遷移確率 0なら(1-c)、1ならcを掛ける
					Tc = 1.0
					for j in range(N*2):
						Tc *= c if ((t >> j) & 1) == 1 else 1.0 - c
					
					dp[i+1][h] = max(dp[i+1][h], dp[i][prev_h] * Tc * E_all)
		
		def trace_back(dps: list[list[float]]) -> list[int]:
			hs = [ 0 ] * M
			h, p = max(enumerate(dp[-1]), key=lambda v: v[1])
			hs[-1] = h
			hp = h >> (N*2)
			hc = h & ((1 << (N*2)) - 1)
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			op = gt_to_int(self.records[-1].v[non_phased_col])
			Ep = E[hp][op]			# parent emission
			for i in range(M-1, 0, -1):
				# どこから遷移してきたのか
				record = self.records[i]
				c = C[i]
				phased_parent_gt = phased_gt_to_int(record.v[phased_col])
				ocs = [ gt_to_int(gt) for gt in record.v[11:] ]
				op = gt_to_int(record.v[non_phased_col])
				E_all = emission_probability(h, op, ocs, phased_parent_gt)
				
				min_diff = 1.0
				opt_h_prev = 0
				for h_prev in range(L):
					t = h_prev ^ hc		# 乗り換えしていないかしたかをbitで
					# 遷移確率 0なら(1-c)、1ならcを掛ける
					Tc = 1.0
					for j in range(N*2):
						Tc *= c if ((t >> j) & 1) == 1 else 1.0 - c
					
					diff = abs(p - dp[i][h_prev] * Tc * E_all)
					if diff < min_diff:
						min_diff = diff
						opt_h_prev = h_prev
				h = opt_h_prev
				hp = h >> (N*2)
				hc = h & ((1 << (N*2)) - 1)
				hs[i-1] = h
				p = dp[i][h]
			
			return hs
		
		def update_genotypes(hs: list[int]):
			phased_col = 9 if self.is_mat_imputed else 10
			non_phased_col = 10 if self.is_mat_imputed else 9
			for i in range(M):
				record = self.records[i]
				hp = hs[i] >> (N*2)
				record.v[non_phased_col] = int_to_phased_gt(hp)
				phased_parent_gt = phased_gt_to_int(record.v[phased_col])
				for j in range(N):	# 個々の後代
					hc = (hs[i] >> (j*2)) & 3
					gtc_int = gt_by_haplotype(hc, hp, phased_parent_gt)
					record.v[j+11] = int_to_phased_gt(gtc_int)
		
		N = self.num_progenies()
		M = len(self)		# マーカー数
		L = 1 << (2*N+2)	# 親とハプロタイプの状態数
		
		# crossover values
		C = [ Map.Kosambi((self.cM(i+1) - self.cM(i)) / 100)
													for i in range(M) ]
		
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E = [ [ 1.0-3*w, w, w, w ],
			  [ w, 1.0-3*w, w, w ],
			  [ w, 1.0-3*w, w, w ],
			  [ w, w, 1.0-3*w, w ] ]
		
		# DP
		# 外側から、マーカー、ハプロタイプの状態
		dp = initialize_dp(L, M)
		for i in range(len(self)):
			update_dp(i, dp)
		
		hs = trace_back(dp)
		update_genotypes(hs)
