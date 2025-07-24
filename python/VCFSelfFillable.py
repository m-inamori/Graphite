from __future__ import annotations

# coding: utf-8
# VCFFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCF import *
from VCFSelfHetero import VCFSelfHetero
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from VCFSelfFillableRecord import VCFSelfFillableRecord
from VCFHeteroHomo import *
from SelfGroups import SelfGroups
from SelfRecordSet import SelfRecordSet
from Genotype import Genotype
from option import *
from common import *


#################### VCFSelfFillable ####################

class VCFSelfFillable(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]],
							records: list[VCFSelfFillableRecord]) -> None:
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records: list[VCFSelfFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def modify(self) -> None:
		# FillTypeでrecordを分ける
		groups = SelfGroups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.determine_parents_phasing()
	
	def find_prev_same_type_record(self, i: int, c: int
										) -> Optional[VCFSelfFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i - 1, -1, -1):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	def find_next_same_type_record(self, i: int, c: int
										) -> Optional[VCFSelfFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i + 1, len(self.records)):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_mat(self, i: int) -> None:
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.is_NA(c-9):
				self.impute_NA_mat_each(i, c)
	
	def create_recordset(self, record: Optional[VCFSelfFillableRecord],
								prev_record: Optional[VCFSelfFillableRecord],
								next_record: Optional[VCFSelfFillableRecord]
															) -> SelfRecordSet:
		return SelfRecordSet(record, prev_record, next_record)
	
	def impute_NA_mat_each(self, i: int, c: int) -> None:
		def select_mat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][0]
#			elif prev_record is None or next_record is None:
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_mat_from
			elif prev_record is None:
				return next_mat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_mat_from
			else:
				return next_mat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		recordset = self.create_recordset(record, prev_record, next_record)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_mat_from = recordset.from_which_chrom_prev_mat(prev_gt)
		next_mat_from = recordset.from_which_chrom_next_mat(next_gt)
		mat_froms = unique_list(prev_mat_from, next_mat_from)
		pairs = [ (x, 1) for x in mat_froms if x != 0 ]
		if not pairs:	# 両側ともにmatが無い（まず無い）
			return
		
		mat_from = select_mat(pairs)
		if mat_from == 0:
			return
		
		record.v[c] = record.gt_from_parent(mat_from, 1) + record.v[c][3:]
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_pat(self, i: int) -> None:
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.v[c][:3] == './.':
				self.impute_NA_pat_each(i, c)
	
	def impute_NA_pat_each(self, i: int, c: int) -> None:
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
		recordset = self.create_recordset(record, prev_record, next_record)
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
	
	@staticmethod
	def fill(vcfs: list[VCFSelfHetero],
					records: list[VCFImpSelfRecord]) -> VCFSelfFillable:
		merged_records = VCFSelfFillable.merge_records(vcfs, records)
		vcf = VCFSelfFillable(vcfs[0].header, merged_records)
		vcf.modify()
		return vcf
	
	@staticmethod
	def merge_records(vcfs: list[VCFSelfHetero],
						records: list[VCFImpSelfRecord],
						) -> list[VCFSelfFillableRecord]:
		all_records = [ VCFSelfFillableRecord.convert(record) for record in
						chain((r for vcf in vcfs for r in vcf.records), records)
		]
		all_records.sort(key=lambda record: record.index)
		return all_records
