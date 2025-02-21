from __future__ import annotations

# coding: utf-8
# VCFOneParentImputed.py
# 片親だけimputeされていてもう片親はknown

from VCFFamily import *
from ParentProgenyImputer import *
from Map import *


#################### VCFOneParentImputed ####################

class VCFOneParentImputed(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records: list[VCFFamilyRecord] = records
		self.imputer = ParentProgenyImputer(records, ref_haps,
													is_mat_imputed, map_, 0.01)
		self.is_mat_imputed = is_mat_imputed
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def impute(self) -> None:
		self.imputer.impute()
