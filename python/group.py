# coding: utf-8
# group.py
# 同じFillTypeのRecordの集まり

from __future__ import annotations
from typing import Optional, Iterator, List, Tuple
from itertools import *

from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from VCFFillableRecord import VCFFillableRecord
from RecordSet import RecordSet, RecordSetSmall

Group = Tuple[FillType, List[VCFFillableRecord]]


#################### Groups ####################

class Groups:
	def __init__(self, groups: list[Group]):
		self.groups: list[Group] = groups
	
	def __len__(self) -> int:
		return len(self.groups)
	
	def find_prev_record(self,
						i: int, g: FillType) -> Optional[VCFFillableRecord]:
		key, records = self.groups[i]
		for j in range(i-1, -1, -1):
			key, records = self.groups[j]
			if key == g:
				return records[-1]
		else:
			return None
	
	def find_next_record(self,
						i: int, g: FillType) -> Optional[VCFFillableRecord]:
		key, records = self.groups[i]
		for j in range(i+1, len(self)):
			key, records = self.groups[j]
			if key == g:
				return records[0]
		else:
			return None
	
	def generate_record_sets(self) -> Iterator[RecordSet]:
		for i, (key, records) in enumerate(self.groups):
			if key == FillType.MAT or key == FillType.PAT:
				continue
			prev_mat_record = self.find_prev_record(i, FillType.MAT)
			next_mat_record = self.find_next_record(i, FillType.MAT)
			prev_pat_record = self.find_prev_record(i, FillType.PAT)
			next_pat_record = self.find_next_record(i, FillType.PAT)
			for record in records:
				yield RecordSet(record, prev_mat_record,
							next_mat_record, prev_pat_record, next_pat_record)
	
	def generate_record_sets_small(self) -> Iterator[RecordSet]:
		for i, (key, records) in enumerate(self.groups):
			if key == FillType.MAT or key == FillType.PAT:
				continue
			prev_mat_record = self.find_prev_record(i, FillType.MAT)
			next_mat_record = self.find_next_record(i, FillType.MAT)
			prev_pat_record = self.find_prev_record(i, FillType.PAT)
			next_pat_record = self.find_next_record(i, FillType.PAT)
			for record in records:
				yield RecordSetSmall(record, prev_mat_record,
							next_mat_record, prev_pat_record, next_pat_record)
	
	@staticmethod
	def create(records: list[VCFFillableRecord]) -> Groups:
		groups = []
		for g, v in groupby(records, key=lambda r: r.type):
			groups.append((g, list(v)))
		return Groups(groups)
