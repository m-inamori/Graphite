from __future__ import annotations

# coding: utf-8
# VCFSelfParentImputed.py
# 自殖で親がimputeされている

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from SelfProgenyImputer import SelfProgenyImputer
from Map import *


#################### VCFSelfParentImputed ####################

class VCFSelfParentImputed(VCFGenoBase):
	def __init__(self, samples: list[str], records: list[GenoRecord],
											map_: Map, vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
		self.records = records
		self.imputers = [ SelfProgenyImputer(records, i, map_, 0.01)
										for i in range(self.num_progenies()) ]
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.samples) - 1
	
	def impute(self) -> None:
		for imputer in self.imputers:
			imputer.impute()
