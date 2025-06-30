from __future__ import annotations

# coding: utf-8
# VCFSelfNoImputedRough.py
# 自殖でimputedなサンプルが一つもない

from VCFFamily import *
from SelfParentImputer import *
from SelfProgenyImputer import *
from Map import *


#################### VCFSelfNoImputedRough ####################

class VCFSelfNoImputedRough(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFRecord],
						ref_haps: list[list[int]], map_: Map) -> None:
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records: list[VCFRecord] = records
		self.parent_imputer = SelfParentImputer(records, ref_haps,
															ic, map_, 0.01)
		self.imputers = [ SelfProgenyImputer(records, i, map_, 0.01)
										for i in range(self.num_progenies()) ]
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 1
	
	def impute(self) -> None:
		self.parent_imputer.impute()
		for imputer in self.imputers:
			imputer.impute()
