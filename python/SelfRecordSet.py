# coding: utf-8
# SelfRecordSet.py
# 自分自身と前後のヘテロ親のセット

from __future__ import annotations
from typing import Optional, Iterator
from itertools import product
from math import log

from VCF import VCFRecord
from VCFSelfFillableRecord import VCFSelfFillableRecord
from TypeDeterminer import ParentComb
from Genotype import Genotype
from log import modified_log
from common import unique_list


#################### SelfRecordSet ####################

class SelfRecordSet:
	def __init__(self, r: Optional[VCFSelfFillableRecord],
					   p: Optional[VCFSelfFillableRecord],
					   n: Optional[VCFSelfFillableRecord]) -> None:
		self.record: Optional[VCFSelfFillableRecord] = r
		self.prev_record: Optional[VCFSelfFillableRecord] = p
		self.next_record: Optional[VCFSelfFillableRecord] = n
	
	def records(self) -> list[Optional[VCFSelfFillableRecord]]:
		return [self.record, self.prev_record, self.next_record]
	
	def prev_mat_from(self, i: int) -> int:
		if self.prev_record is None:
			return 0
		return self.prev_record.from_which_chrom(i, True)
	
	def next_mat_from(self, i: int) -> int:
		if self.next_record is None:
			return 0
		return self.next_record.from_which_chrom(i, True)
	
	def prev_pat_from(self, i: int) -> int:
		if self.prev_record is None:
			return 0
		return self.prev_record.from_which_chrom(i, False)
	
	def next_pat_from(self, i: int) -> int:
		if self.next_record is None:
			return 0
		return self.next_record.from_which_chrom(i, False)
	
	def gt_each(self, i: int, record: Optional[VCFSelfFillableRecord]) -> str:
		return './.' if record is None else record.v[i+9]
	
	def gt(self, i: int) -> str:
		return self.gt_each(i, self.record)
	
	def gts(self, i: int) -> tuple[str, str, str]:
		v = [ self.gt_each(i, r) for r in self.records() ]
		return (v[0], v[1], v[2])
	
	def __select_phasing(self, candidates: list[int]) -> int:
		if len(candidates) == 1 or self.record is None:
			return candidates[0]
		
		parent_gt = self.record.get_int_gt(0)
		
		v: list[tuple[float, int]] = []
		for phasing in candidates:
			parent_int_gt1 = (phasing >> 1) + (phasing & 1)
			score = abs(parent_int_gt1 - parent_gt)
			v.append((score, phasing))
		
		# これだと、同じスコアならphasingが小さい方から選んでいる
		# 本当はランダム的に選びたい
		_, phasing = min(v)
		return phasing
	
	def determine_phasing_core(self, lls: list[tuple[float, int]]) -> int:
		# 最大に近いllを集める
		candidates = []
		max_ll = lls[-1][0]
		for ll, phasing in reversed(lls):
			if ll > max_ll - 1e-9:
				candidates.append(phasing)
			else:
				break
		return self.__select_phasing(candidates)
	
	def possible_phasings(self) -> list[int]:
		if self.record is None:
			return []
		
		if self.record.pair == ParentComb.P00x00:
			return [0]
		elif self.record.pair == ParentComb.P11x11:
			return [3]
		else:
			# all
			return [0, 1, 2, 3]
	
	# phasingされている前提
	def from_which_chrom(self, gt: str, record: Optional[VCFRecord],
														mat: bool) -> int:
		if record is None:
			return 0
		
		i = 0 if mat else 1
		parent_gt = record.v[9]
		return 1 if parent_gt[0] == gt[i*2] else 2
	
	def from_which_chrom_prev_mat(self, gt: str) -> int:
		return self.from_which_chrom(gt, self.prev_record, True)
	
	def from_which_chrom_next_mat(self, gt: str) -> int:
		return self.from_which_chrom(gt, self.next_record, True)
	
	def from_which_chrom_prev_pat(self, gt: str) -> int:
		return self.from_which_chrom(gt, self.prev_record, False)
	
	def from_which_chrom_next_pat(self, gt: str) -> int:
		return self.from_which_chrom(gt, self.next_record, False)
	
	def gen_gts(self) -> Iterator[tuple[int, str, str]]:
		if self.record is None:
			return
		
		for c in range(11, len(self.record.v)):
			gts = [ r.v[c] if r else '' for r in self.records() ]
			yield (c-9, gts[1], gts[2])
	
	def is_prev_near(self) -> bool:
		if (self.record is None or self.prev_record is None or
										self.next_record is None):
			return False
		return (self.record.pos() * 2 <
					self.prev_record.pos() + self.next_record.pos())
	
	def near_mat_from(self, i: int) -> int:
		return (self.prev_mat_from(i) if self.is_prev_near()
										else self.next_mat_from(i))
	
	def near_pat_from(self, i: int) -> int:
		return (self.prev_pat_from(i) if self.is_prev_near()
										else self.next_pat_from(i))
	
	# [(mat_from, pat_from)] -> (mat_from, pat_from)
	def select_nearest_froms(self, pairs: list[tuple[int, int]],
											i: int) -> tuple[int, int]:
		if len(pairs) == 4:		# N/Aのときにまれにあり得る
			return (self.near_mat_from(i), self.near_pat_from(i))
		elif pairs[0][0] == pairs[1][0]:		# matが同じ
			if (self.record is None or self.prev_record is None or
											self.next_record is None):
				return (0, 0)
			elif self.is_prev_near():
				return (pairs[0][0], self.prev_pat_from(i))
			else:
				return (pairs[0][0], self.next_pat_from(i))
		elif pairs[0][1] == pairs[1][1]:	# patが同じ
			if (self.record is None or self.prev_record is None or
											self.next_record is None):
				return (0, 0)
			elif self.is_prev_near():
				return (self.prev_mat_from(i), pairs[0][1])
			else:
				return (self.next_mat_from(i), pairs[0][1])
		else:	# 両親とも乗り換えている（滅多にない）
			return (self.near_mat_from(i), self.near_pat_from(i))
	
	def select_pair(self, pairs: list[tuple[int, int]], i: int,
							selected: bool = False) -> tuple[int, int]:
		def sum_gt(gt: str) -> int:
			try:
				return int(gt[0]) + int(gt[2])
			except ValueError:
				return -1
		
		record = self.record
		if record is None:	# for mypy
			return (0, 0)
		parent_gt = record.get_int_gt(0)
		gt = record.get_gt(i)
		if not pairs:
			return (0, 0)
		elif len(pairs) == 1:
			return pairs[0]
		elif not Genotype.is_valid(gt, parent_gt, parent_gt):
			return self.select_nearest_froms(pairs, i)
		elif selected:
			return self.select_nearest_froms(pairs, i)
		else:
			new_pairs = [ v for v in pairs
						if sum_gt(gt) == sum_gt(record.gt_from_parent(*v)) ]
			pair = self.select_pair(new_pairs, i, True)
			if pair != (0, 0):
				return pair
			else:
				return self.select_pair(pairs, i, True)
	
	# mat_gt, pat_gt : 0|0 0|1 1|0 1|1を0～3で表す
	def compute_phasing_likelihood(self, phasing: int) -> float:
		memo = { (0, 0): [0.5, 0.5], (0, 1): [0.9, 0.1], (0, 2): [0.1, 0.9],
				 (1, 0): [0.9, 0.1], (1, 1): [0.99, 0.01], (1, 2): [0.5, 0.5],
				 (2, 0): [0.1, 0.9], (2, 1): [0.5, 0.5], (2, 2): [0.01, 0.99] }
		def probs_from_which_chrom(prev_chrom: int,
								   next_chrom: int) -> list[float]:
			return memo[(prev_chrom, next_chrom)]
		
		def likelihood_each(probs_mat: list[float],
							probs_pat: list[float], i: int) -> float:
			if self.record is None:
				return log(0.0001)
			
			likelihood = 0.0
			for j, k in product(range(2), repeat=2):
				gt = ((phasing >> j) & 1) + ((phasing >> k) & 1)
				likelihood += (probs_mat[j] * probs_pat[k] *
												self.record.probs[i][gt])
			return modified_log(likelihood)
		
		if self.record is None:
			return log(0.0001)
		
		ll = 0.0
		for i, gt1, gt2 in self.gen_gts():
			prev_mat_from = self.from_which_chrom_prev_mat(gt1)
			next_mat_from = self.from_which_chrom_next_mat(gt2)
			prev_pat_from = self.from_which_chrom_prev_pat(gt1)
			next_pat_from = self.from_which_chrom_next_pat(gt2)
			probs_mat = probs_from_which_chrom(prev_mat_from, next_mat_from)
			probs_pat = probs_from_which_chrom(prev_pat_from, next_pat_from)
			ll += likelihood_each(probs_mat, probs_pat, i)
		return ll
	
	def determine_parents_phasing(self) -> None:
		if self.record is None:
			return
		
		record = self.record
		lls = []
		for phasing in self.possible_phasings():
			ll = self.compute_phasing_likelihood(phasing)
			lls.append((ll, phasing))
		lls.sort()
		
		phasing = self.determine_phasing_core(lls)
		if record is None:
			return
		gt = ['0|0', '1|0', '0|1', '1|1']
		record.set_GT(0, gt[phasing])
