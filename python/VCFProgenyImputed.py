from __future__ import annotations

# coding: utf-8
# VCFProgenyImputed.py
# 片親がknownで後代がimputeされている家系を補完する
# 2サンプルのみ

from VCFFamily import *
from ParentImputerByProgeny import *


#################### VCFProgenyImputed ####################

class VCFProgenyImputed(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFRecord],
					ref_haps: list[list[int]], is_mat_known: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.imputer = ParentImputerByProgeny(records, ref_haps,
												is_mat_known, map_, 0.01)
		self.records: list[VCFRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2
	
	def impute(self) -> None:
		self.imputer.impute()
