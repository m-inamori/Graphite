from __future__ import annotations

# coding: utf-8
# VCFOneParentImputedRough.py
# 片側だけimputeされている
# 親だけでDPで補完する
# 子でもペナルティをつける
# 例えば、親のGenotypeが0|1で子に0|0があったらペナルティ


from VCFFamily import *
from ParentImputer import *
from ProgenyImputer import *
from Map import *


#################### VCFOneParentImputedRough ####################

class VCFOneParentImputedRough(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records = records
		self.parent_imputer = ParentImputer(records, not is_mat_imputed,
														ref_haps, map_, 0.01)
		self.prog_imputer = ProgenyImputer(records, map_, 0.01)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2
	
	def impute(self) -> None:
		self.parent_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
