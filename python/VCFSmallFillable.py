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
