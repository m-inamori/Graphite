from __future__ import annotations

# coding: utf-8
# VCFOneParentImputed.py
# 片親だけimputeされていてもう片親はknown

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFFamily import *
from ParentProgenyImputer import *
from Map import *


#################### VCFOneParentImputed ####################

class VCFOneParentImputed(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool,
					map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.imputer = ParentProgenyImputer(records, ref_haps,
													is_mat_imputed, map_, 0.01)
		self.is_mat_imputed = is_mat_imputed
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> None:
		self.imputer.impute()
