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
