from __future__ import annotations

# coding: utf-8
# VCFSelfNoImputed.py
# 自殖でimputedなサンプルが一つもない

from VCFFamily import *
from SelfParentImputer import *
from SelfProgenyImputer import *
from Map import *


#################### VCFSelfNoImputed ####################

class VCFSelfNoImputed(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFRecord],
						ref_haps: list[list[int]], map_: Map) -> None:
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records: list[VCFRecord] = records
		self.imputer = SelfImputer(records, ref_haps, map_, 0.01)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 1
	
	def impute(self) -> None:
		self.imputer.impute()
