from __future__ import annotations

# coding: utf-8
# VCFHeteroImpHomo.py
# phasingされている方がホモのヘテロ×ホモファミリー

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from VCFFillable import *
from group import Groups
from RecordSet import RecordSet
from VCFImpFamily import FillType
from VCFHeteroHomoPP import *
from VCFHeteroHomoOnePhased import *
from Map import *
import ClassifyRecord
import Imputer
from option import *


#################### State ####################

class State:
	def __init__(self, s: int, n: int):
		self.s = s
		self.n = n
	
	def order(self) -> int:
		return self.s & 1
	
	def haplotype(self, i: int) -> int:
		return (self.s >> (i+1)) & 1
	
	def num_crossovers(self) -> int:
		return self.s >> (self.n+1)
	
	def set_haplotype(self, h: int):
		self.s |= h << (self.n+1)
	
	def increment_num_crossovers(self):
		self.s += 1 << (self.n+1)


#################### VCFHeteroImpHomo ####################

class VCFHeteroImpHomo(VCFHeteroHomoOnePhased):
	def __init__(self, header: list[list[str]],
						records: list[VCFFillableRecord],
						is_mat_hetero: bool, map_: Map):
		VCFHeteroHomoOnePhased.__init__(self, header, records,
												is_mat_hetero, map_)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def imputed_index(self):
		return 1 if self.is_mat_hetero else 0
	
	def unimputed_index(self):
		return 0 if self.is_mat_hetero else 1
	
	def impute(self):
		INF = 10**9
		def init_dp(n: int) -> list[tuple[int, int]]:
			dp = [(INF, 0)] * ((n*2+1)*2**(n+1))
			dp[0] = (0, 0)
			return dp
		
		# (親は0|1と1|0のどちらか, 総乗り換え数, 
		def decode(stat: int, n: int) -> tuple[int, list[int], int]:
			v = []
			stat, order = divmod(stat, 2)
			for _ in range(n):
				stat, r = divmod(stat, 2)
				v.append(r)
			return (order, v, stat)
		
		# ハプロタイプを決めたときのGenotype
		def gt_by_haplotype(h: int, gt_imputed: int, order: int) -> int:
			return gt_imputed // 2 + (order ^ h)
		
		# ヘテロ親が決まっている時、各後代が乗り換えしないでも正しいか
		def is_right_gt(order: int, gt_imputed: int, prev_hs: list[int],
											record: VCFRecord) -> list[bool]:
			bs: list[bool] = []
			for i in range(2, n+2):
				gt = record.get_int_gt(i)
				if gt == -1:
					# N/Aならミス無しとみなして乗り換えしない
					bs.append(True)
				else:
					new_gt = gt_by_haplotype(prev_hs[i-2], gt_imputed, order)
					bs.append(new_gt == gt)
			return bs
		
		def next_states(state: State,
							record: VCFRecord) -> Iterator[tuple[State, int]]:
			gt_imputed = record.get_int_gt(self.imputed_index())
			# order 0 -> '0|1', 1 -> '1|0'
			for order in (0, 1):
				bs = is_right_gt(order, gt_imputed, prev_hs, record)
				num_falses = sum(1 for b in bs if not b)
				s = State(order, n)
				s.set_num_crossovers(state.num_crossovers())
				mis = 0
				for s1 in range(1 << num_falses):
					j = 0
					for i in range(2, n+2):
						prev_h = prev_hs[i-2]
						if not bs[i-2]:		# 乗り換えなしだと違う遺伝子型になる
							h = (s1 >> j) & 1
							if h != prev_h:		# 乗り換えあり
								if s >> (n+1) < n * 2:
									s |= h << (i-1)
									s += 1 << (n+1)
							else:				# 乗り換えなしだが違う遺伝子型
								mis += 1
							j += 1
						else:
							h = prev_h
							s |= h << (i-1)
					yield (s, mis)
		
		def update_dp(dp: list[tuple[int, int]],
							record: VCFRecord) -> list[tuple[int, int]]:
			new_dp = [(INF, 0)] * len(dp)
			n = len(record.samples) - 2
			for state, (mis, _) in enumerate(dp):
				for new_state, new_mis in next_states(state, record):
					value = (mis+new_mis, state)
					new_dp[new_state] = min(new_dp[new_state], value)
			return new_dp
		
		def trace_back(state: int, dps: list[list[tuple[int, int]]]):
			n = self.num_progenies()
			for i in range(len(self)-1, -1, -1):
				order, hs, _ = decode(state, n)
				record = self.records[i]
				hetero_GT = '0|1' if order == 0 else '1|0'
				record.set_GT(self.unimputed_index(), hetero_GT)
				gt_imp = record.get_gt(self.imputed_index())[0]
				for j in range(2, len(self.samples)):
					h = hs[j-2]
					if self.is_mat_hetero:
						GT = "%d|%s" % (order ^ h, gt_imp)
					else:
						GT = "%s|%d" % (gt_imp, order ^ h)
					record.set_GT(j, GT)
				_, state = dps[i+1][state]
		
		n = self.num_progenies()
		dps = [init_dp(n)]
		for record in self.records:
			dps.append(update_dp(dps[-1], record))
		_, state = min(zip(dps[-1], count()))
		trace_back(state, dps)
