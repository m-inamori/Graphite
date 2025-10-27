from __future__ import annotations

# coding: utf-8
# VCFProgenyImputed.py
# 片親がknownで後代がimputeされている家系を補完する
# 2サンプルのみ

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from ParentImputerByProgeny import *


#################### VCFProgenyImputed ####################

class VCFProgenyImputed(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]],
					is_mat_known: bool, map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		self.imputer = ParentImputerByProgeny(records, ref_haps,
												is_mat_known, map_, 0.01)
		self.records: list[VCFFamilyRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> None:
		self.imputer.impute()
