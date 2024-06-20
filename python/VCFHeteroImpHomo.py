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
		self.s = s	# haplotypes & num_crossovers
		self.n = n		# num of pogenies
	
	def haplotype(self, i: int) -> int:
		return (self.s >> i) & 1
	
	def num_crossovers(self) -> int:
		return self.s >> self.n
	
	def is_full_num_crossovers(self) -> bool:
		return self.num_crossovers() >= self.n * 2
	
	def set_haplotype(self, h: int, i: int):
		self.s |= h << i
	
	def increment_num_crossovers(self):
		self.s += 1 << self.n
	
	def set_num_crossovers(self, n):
		self.s |= n << self.n
	
	@staticmethod
	def create_default(n: int) -> State:
		return State(0, n)


#################### Value ####################

INF = 10**9

class Value:
	def __init__(self, mis: int, state: State, o: int):
		self.num_mis: int = mis
		self.state: State = state
		self.order: int = o
	
	def is_valid(self) -> bool:
		return self.num_mis != INF
	
	def add(self, other: Value, state: State) -> Value:
		return Value(self.num_mis + other.num_mis, state, other.order)
	
	def __lt__(self, other: Value) -> bool:
		if self.num_mis != other.num_mis:
			return self.num_mis < other.num_mis
		elif self.state.num_crossovers() != other.state.num_crossovers():
			return self.state.num_crossovers() < other.state.num_crossovers()
		else:
			return self.order < other.order
	
	@staticmethod
	def create_default(n: int) -> Value:
		return Value(INF, State.create_default(n), 0)


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
		def init_dp(n: int) -> list[Value]:
			dp = [ Value.create_default(n) for _ in range((n*2+1)*2**n) ]
			dp[0] = Value(0, State(0, n), 0)
			return dp
		
		# ハプロタイプを決めたときのGenotype
		def gt_by_haplotype(h: int, gt_imputed: int, order: int) -> int:
			return gt_imputed // 2 + (order ^ h)
		
		# ヘテロ親が決まっている時、各後代が乗り換えしないでも正しいか
		def is_right_gt(order: int, gt_imputed: int, state: State,
											record: VCFRecord) -> list[bool]:
			bs: list[bool] = []
			for i in range(n):
				gt = record.get_int_gt(i+2)
				if gt == -1:
					# N/Aならミス無しとみなして乗り換えしない
					bs.append(True)
				else:
					prev_h = state.haplotype(i)
					new_gt = gt_by_haplotype(prev_h, gt_imputed, order)
					bs.append(new_gt == gt)
			return bs
		
		def next_states(state: State,
						record: VCFRecord) -> Iterator[tuple[State, Value]]:
			n = self.num_progenies()
			gt_imputed = record.get_int_gt(self.imputed_index())
			# order 0 -> '0|1', 1 -> '1|0'
			for order in (0, 1):
				bs = is_right_gt(order, gt_imputed, state, record)
				# 乗り換えなしで遺伝子型が正しいか
				num_falses = sum(1 for b in bs if not b)
				s0 = State(0, n)
				s0.set_num_crossovers(state.num_crossovers())
				for s1 in range(1 << num_falses):
					s = State(s0.s, n)
					mis = 0
					j = 0	# num_falsesの分だけインクリメントする
					for i in range(n):
						prev_h = state.haplotype(i)
						if not bs[i]:		# 乗り換えなしだと違う遺伝子型になる
							h = (s1 >> j) & 1
							if h != prev_h:		# 乗り換えあり
								if s.is_full_num_crossovers():
									break
								s.increment_num_crossovers()
							else:				# 乗り換えなしだが違う遺伝子型
								mis += 1
							j += 1
						else:
							h = prev_h	# s1と無関係に乗り換えない
						s.set_haplotype(h, i)
					else:
						yield (s, Value(mis, state, order))
		
		def update_dp(dp: list[Value], record: VCFRecord) -> list[Value]:
			n = len(record.samples) - 2
			new_dp = [ Value.create_default(n) for _ in range(len(dp)) ]
			for s, value0 in enumerate(dp):
				if not value0.is_valid():
					continue
				state = State(s, n)
				for new_state, new_value in next_states(state, record):
					value = value0.add(new_value, state)
					new_dp[new_state.s] = min(new_dp[new_state.s], value)
			return new_dp
		
		def trace_back(state: State, dps: list[list[Value]]):
			n = self.num_progenies()
			for i in range(len(self)-1, -1, -1):
				order = dps[i+1][state.s].order
				record = self.records[i]
				hetero_GT = '0|1' if order == 0 else '1|0'
				record.set_GT(self.unimputed_index(), hetero_GT)
				gt_imp = record.get_gt(self.imputed_index())[0]
				for j in range(2, len(self.samples)):
					h = state.haplotype(j-2)
					if self.is_mat_hetero:
						GT = "%d|%s" % (order ^ h, gt_imp)
					else:
						GT = "%s|%d" % (gt_imp, order ^ h)
					record.set_GT(j, GT)
				state = dps[i+1][state.s].state
		
		n = self.num_progenies()
		dps = [init_dp(n)]
		for record in self.records:
			dps.append(update_dp(dps[-1], record))
		
		_, s = min(zip(dps[-1], count()), key=lambda v: (v[0].num_mis, v[1]))
		state = State(s, n)
		trace_back(state, dps)
