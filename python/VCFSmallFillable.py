from __future__ import annotations

# coding: utf-8
# VCFSmallFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFGeno import VCFGeno
from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
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
	def __init__(self, samples: list[str],
							records: list[VCFFillableRecord], vcf: VCFSmall):
		VCFFillable.__init__(self, samples, records, vcf)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
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
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_pat(self, i: int) -> None:
		record = self.records[i]
		for k in range(2, len(record.geno)):
			if record.is_NA(k):
				self.impute_NA_pat_each(i, k)
	
	def remove_parents(self) -> VCFGeno:
		samples: list[str] = self.samples[2:]
		records: list[GenoRecord] = []
		for record in self.records:
			geno = record.geno[2:]
			r = GenoRecord(record.pos, geno)
			records.append(r)
		return VCFGeno(samples, records, self.vcf)


__all__ = ['VCFSmallFillable']
