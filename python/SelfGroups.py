# coding: utf-8
# SelfGroups.py
# 同じSelfFillTypeのRecordの集まり

from __future__ import annotations
from typing import Optional, Iterator, List, Tuple
from itertools import *

from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from VCFSelfFillableRecord import VCFSelfFillableRecord
from SelfRecordSet import SelfRecordSet

Group = Tuple[SelfFillType, List[VCFSelfFillableRecord]]


#################### SelfGroups ####################

class SelfGroups:
	def __init__(self, groups: list[Group]):
		self.groups: list[Group] = groups
	
	def __len__(self) -> int:
		return len(self.groups)
	
	def find_prev_record(self, i: int, g: SelfFillType
								) -> Optional[VCFSelfFillableRecord]:
		key, records = self.groups[i]
		chr = records[0].v[0]
		for j in range(i-1, -1, -1):
			key, records = self.groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[-1]
		else:
			return None
	
	def find_next_record(self, i: int, g: SelfFillType
									) -> Optional[VCFSelfFillableRecord]:
		key, records = self.groups[i]
		chr = records[0].v[0]
		for j in range(i+1, len(self)):
			key, records = self.groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[0]
		else:
			return None
	
	def generate_record_sets(self) -> Iterator[SelfRecordSet]:
		for i, (key, records) in enumerate(self.groups):
			if key == SelfFillType.P01 or key == SelfFillType.FILLED:
				continue
			prev_record = self.find_prev_record(i, SelfFillType.P01)
			next_record = self.find_next_record(i, SelfFillType.P01)
			for record in records:
				yield SelfRecordSet(record, prev_record, next_record)
	
	@staticmethod
	def create(records: list[VCFSelfFillableRecord]) -> SelfGroups:
		groups: list[Group] = []
		for g, v in groupby(records, key=lambda r: r.type):
			groups.append((g, list(v)))
		return SelfGroups(groups)
