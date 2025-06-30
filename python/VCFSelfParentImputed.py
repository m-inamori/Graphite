from __future__ import annotations

# coding: utf-8
# VCFSelfParentImputed.py
# 自殖で親がimputeされている

from VCFFamily import *
from SelfProgenyImputer import *
from Map import *


#################### VCFSelfParentImputed ####################

class VCFSelfParentImputed(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFRecord], map_: Map) -> None:
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records = records
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
		for imputer in self.imputers:
			imputer.impute()
