from __future__ import annotations

# coding: utf-8
# VCFSmallFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from VCFFillableRecord import VCFFillableRecord
from VCFFillable import *
from VCFHeteroHomo import *
from group import Groups
from RecordSet import RecordSet, RecordSetSmall
from Genotype import Genotype
from option import *
from common import *


#################### VCFSmallFillable ####################

class VCFSmallFillable(VCFFillable):
	def __init__(self, header: list[list[str]],
							records: list[VCFFillableRecord]):
		VCFFillable.__init__(self, header, records)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def modify(self, is_phased_changable: bool) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets_small():
			record_set.determine_parents_phasing()
			self.impute_core(record_set)
		
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == FillType.MAT:
				self.impute_NA_mat(i)
			elif record.type == FillType.PAT:
				self.impute_NA_pat(i)
			elif record.type in (FillType.IMPUTABLE, FillType.UNABLE):
				self.impute_others(i)
		
		for record in self.records:
			record.fill_PGT()
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_pat(self, i: int):
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.is_NA(c-9):
				self.impute_NA_pat_each(i, c)
	
	"""
	def impute_NA_pat_each(self, i: int, c: int):
		def select_pat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][1]
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_pat_from
			elif prev_record is None:
				return next_pat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_pat_from
			else:
				return next_pat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		recordset = RecordSetSmall(record, None, None, prev_record, next_record)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_pat_from = recordset.from_which_chrom_prev_pat(prev_gt)
		next_pat_from = recordset.from_which_chrom_next_pat(next_gt)
		pat_froms = unique_list(prev_pat_from, next_pat_from)
		pairs = [ (1, x) for x in pat_froms if x != 0 ]
		if not pairs:	# 両側ともにpatが無い（まず無い）
			return
		
		pat_from = select_pat(pairs)
		record.v[c] = record.gt_from_parent(1, pat_from) + record.v[c][3:]
	"""
	
	"""
	def __impute_others(self, i: int):
		def correct(GT: str, c: int, record: VCFFillableRecord) -> str:
			int_gt = record.get_int_gt(c-9)
			if int_gt == -1:
				if GT[0] != '.':
					return GT[:2] + str(record.pos() % 2)
				elif GT[2] != '.':
					return str(record.pos() % 2) + GT[1:]
				else:
					return (str(record.pos() % 2) + '|' +
									str((record.pos() >> 1) % 2))
			elif int_gt == 0:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					return '0|0'
			elif int_gt == 1:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					# どちらにすべきか判別がつかないので適当に選ぶ
					return '0|1' if int(record.v[1]) % 2 == 0 else '1|0'
			else:
				if GT[0] == '0':
					return '0|1'
				elif GT[2] == '0':
					return '1|0'
				else:
					return '1|1'
		
		record = self.records[i]
		cs = [ c for c in range(11, len(record.v))
				if record.v[c][1] == '/' or record.get_int_gt(c-9) == -1 ]
		if not cs:
			return
		
		mat_homo = record.is_homo(0)
		pat_homo = record.is_homo(1)
		for c in cs:
			mat_from = 1 if mat_homo else self.__find_mat_from(i, c)
			pat_from = 1 if pat_homo else self.__find_pat_from(i, c)
			GT = record.gt_from_parent(mat_from, pat_from)
			if '.' in GT:
				GT = correct(GT, c, record)
			record.set_GT(c-9, GT)
	
	def __find_mat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_mat_from(i, c),
									self.__find_next_mat_from(i, c), i)
	
	def __find_pat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_pat_from(i, c),
									self.__find_next_pat_from(i, c), i)
	
	def __select_from(self, f1: tuple[int, int],
							f2: tuple[int, int], i: int) -> int:
		i1, from1 = f1
		i2, from2 = f2
		if from1 == 0 and from2 == 0:
			# 前後がないとき乱数的に決める
			r0 = self.records[i]
			return r0.pos() % 2 + 1
		
		if from1 == from2:
			return from1
		elif from2 == 0:
			return from1
		elif from1 == 0:
			return from2
		else:
			# 最後は物理距離で決める
			r0 = self.records[i]
			r1 = self.records[i1]
			r2 = self.records[i2]
			if r0.pos() * 2 <= r2.pos() + r1.pos():
				return from1
			else:
				return from2
	
	def __find_prev_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_prev_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	"""
	
	def remove_parents(self) -> VCFSmall:
		samples: list[str] = self.samples[2:]
		records: list[VCFRecord] = []
		for record in self.records:
			v = record.v[:9] + record.v[11:]
			r = VCFRecord(v, samples)
			records.append(r)
		header = self.trim_header(samples)
		return VCFSmall(header, records)


__all__ = ['VCFSmallFillable']
