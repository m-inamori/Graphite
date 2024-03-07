# coding: utf-8
# RecordSet.py
# 自分自身と前後のヘテロ×ホモ親のセット

from __future__ import annotations
from typing import Optional, Iterator
from itertools import product
from math import log

from VCF import VCFRecord
from VCFFillableRecord import VCFFillableRecord
from TypeDeterminer import ParentComb
from Genotype import Genotype


#################### RecordSet ####################

class RecordSet:
	def __init__(self, r: Optional[VCFFillableRecord],
			pm: Optional[VCFFillableRecord], nm: Optional[VCFFillableRecord],
			pp: Optional[VCFFillableRecord], np: Optional[VCFFillableRecord]):
		self.record: Optional[VCFFillableRecord] = r
		self.prev_mat_record: Optional[VCFFillableRecord] = pm
		self.next_mat_record: Optional[VCFFillableRecord] = nm
		self.prev_pat_record: Optional[VCFFillableRecord] = pp
		self.next_pat_record: Optional[VCFFillableRecord] = np
	
	def records(self) -> list[Optional[VCFFillableRecord]]:
		return [self.record, self.prev_mat_record, self.next_mat_record,
							 self.prev_pat_record, self.next_pat_record]
	
	def gt_each(self, i: int, record: Optional[VCFFillableRecord]) -> str:
		return './.' if record is None else record.v[i+9]
	
	def gt(self, i: int) -> str:
		return self.gt_each(i, self.record)
	
	def gts(self, i: int) -> tuple[str, str, str, str, str]:
		v = [ self.gt_each(i, r) for r in self.records() ]
		return (v[0], v[1], v[2], v[3], v[4])
	
	def is_prev_nearer(self, is_mat: bool) -> bool:
		if self.record is None:
			return False
		elif is_mat:
			# どちらかがNoneの場合、ここには来ないので、適当に処理する
			if self.prev_mat_record is None or self.next_mat_record is None:
				return False
			else:
				return (self.prev_mat_record.pos() + self.next_mat_record.pos()
														>= self.record.pos()*2)
		else:
			if self.prev_pat_record is None or self.next_pat_record is None:
				return False
			else:
				return (self.prev_pat_record.pos() + self.next_pat_record.pos()
														>= self.record.pos()*2)
	
	def __select_phasing(self, candidates: list[tuple[int, int]]
												) -> tuple[int, int]:
		if len(candidates) == 1 or self.record is None:
			return candidates[0]
		
		mat_int_gt = self.record.mat_int_gt()
		pat_int_gt = self.record.pat_int_gt()
		
		v: list[tuple[float, int, int]] = []
		for mat_phasing, pat_phasing in candidates:
			mat_int_gt1 = (mat_phasing >> 1) + (mat_phasing & 1)
			pat_int_gt1 = (pat_phasing >> 1) + (pat_phasing & 1)
			score = (abs(mat_int_gt1 - mat_int_gt) +
					 abs(pat_int_gt1 - pat_int_gt))
			v.append((score, mat_phasing, pat_phasing))
		
		# これだと、同じスコアならphasingが小さい方から選んでいる
		# 本当はランダム的に選びたい
		_, mat_phasing, pat_phasing = min(v)
		return (mat_phasing, pat_phasing)
	
	def determine_phasing_core(self, lls: list[tuple[float, int, int]]
														) -> tuple[int, int]:
		# 最大に近いllを集める
		candidates = []
		max_ll = lls[-1][0]
		for ll, mat_phasing, pat_phasing in reversed(lls):
			if ll > max_ll - 1e-9:
				candidates.append((mat_phasing, pat_phasing))
			else:
				break
		return self.__select_phasing(candidates)
	
	def possible_phasings(self) -> list[tuple[int, int]]:
		# pair 0: 0x0, 1: 0x1, 2: 0x2, 3: 1x1, 4: 1x2, 5: 2x2
		if self.record is None:
			return []
		elif self.record.pair == ParentComb.P00x00:
			return [(0, 0)]
		elif self.record.pair == ParentComb.P00x01:
			return [(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.record.pair == ParentComb.P01x01:
			return [(1, 1), (1, 2), (2, 1), (2, 2),
					(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 3), (2, 3), (3, 1), (3, 2)]
		elif self.record.pair == ParentComb.P00x11:
			return [(0, 3), (3, 0)]
		elif self.record.pair == ParentComb.P01x11:
			return [(1, 3), (2, 3), (3, 1), (3, 2),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.record.pair == ParentComb.P11x11:
			return [(3, 3)]
		else:
			# all
			return list(product(range(4), range(4)))
	
	# phasingされている前提
	def from_which_chrom(self, gt: str, record: Optional[VCFRecord],
														mat: bool) -> int:
		if record is None:
			return 0
		
		i = 0 if mat else 1
		parent_gt = record.v[9+i]
		return 1 if parent_gt[0] == gt[i*2] else 2
	
	def from_which_chrom_prev_mat(self, gt) -> int:
		return self.from_which_chrom(gt, self.prev_mat_record, True)
	
	def from_which_chrom_next_mat(self, gt) -> int:
		return self.from_which_chrom(gt, self.next_mat_record, True)
	
	def from_which_chrom_prev_pat(self, gt) -> int:
		return self.from_which_chrom(gt, self.prev_pat_record, False)
	
	def from_which_chrom_next_pat(self, gt) -> int:
		return self.from_which_chrom(gt, self.next_pat_record, False)
	
	def gen_gts(self) -> Iterator[tuple[str, str, str, str, str]]:
		if self.record is None:
			return
		for c in range(11, len(self.record.v)):
			gts = [ r.v[c] if r else '' for r in self.records() ]
			yield (gts[0], gts[1], gts[2], gts[3], gts[4])
	
	# mat_gt, pat_gt : 0|0 0|1 1|0 1|1を0～3で表す
	def compute_phasing_likelihood(self, mat_gt: int, pat_gt: int) -> float:
		memo = { (0, 0): [0.5, 0.5], (0, 1): [0.9, 0.1], (0, 2): [0.1, 0.9],
				 (1, 0): [0.9, 0.1], (1, 1): [0.99, 0.01], (1, 2): [0.5, 0.5],
				 (2, 0): [0.1, 0.9], (2, 1): [0.5, 0.5], (2, 2): [0.01, 0.99] }
		def probs_from_which_chrom(prev_chrom: int,
								   next_chrom: int) -> list[float]:
			return memo[(prev_chrom, next_chrom)]
		
		def likelihood_each(gt: str, probs_mat: list[float],
									 probs_pat: list[float]) -> float:
			sum_gt = int(gt[0]) + int(gt[2])
			likelihood = sum(probs_mat[i] * probs_pat[j]
								for i, j in product(range(2), repeat=2)
								if ((mat_gt >> i) & 1) +
								   ((pat_gt >> j) & 1) == sum_gt)
			return log(likelihood)
		
		ll = 0.0	# log of likelihood
		for gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in self.gen_gts():
			if not Genotype.is_valid(gt, mat_gt, pat_gt):
				ll += log(0.0001)
				continue
			prev_mat_from = self.from_which_chrom_prev_mat(mat_gt1)
			next_mat_from = self.from_which_chrom_next_mat(mat_gt2)
			prev_pat_from = self.from_which_chrom_prev_pat(pat_gt1)
			next_pat_from = self.from_which_chrom_next_pat(pat_gt2)
			probs_mat = probs_from_which_chrom(prev_mat_from, next_mat_from)
			probs_pat = probs_from_which_chrom(prev_pat_from, next_pat_from)
			ll += likelihood_each(gt, probs_mat, probs_pat)
		return ll
	
	def determine_parents_phasing(self):
		if self.record is None:
			return
		
		record = self.record
		lls = []
		for mat_p, pat_p in self.possible_phasings():
			ll = self.compute_phasing_likelihood(mat_p, pat_p)
			lls.append((ll, mat_p, pat_p))
		lls.sort()
		
		mat_gt, pat_gt = self.determine_phasing_core(lls)
		if record is None:
			return
		gt = ['0|0', '1|0', '0|1', '1|1']
		record.v[9] = gt[mat_gt] + record.v[9][3:]
		record.v[10] = gt[pat_gt] + record.v[10][3:]
