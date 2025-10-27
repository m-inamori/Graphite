from __future__ import annotations

# coding: utf-8
# VCFSelfProgenyImputed.py
# 自殖で後代の一部がimputeされている

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from SelfParentImputer import SelfParentImputer
from SelfProgenyImputer import SelfProgenyImputer
from Map import *


#################### VCFSelfProgenyImputed ####################

class VCFSelfProgenyImputed(VCFGenoBase, VCFMeasurable):
	def __init__(self, samples: list[str],
						records: list[GenoRecord],
						ref_haps: list[list[int]],
						ic: int, map_: Map, vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records: list[GenoRecord] = records
		self.ic: int = ic	# index of imputed progeny
		self.parent_imputer = SelfParentImputer(records, ref_haps,
															ic, map_, 0.01)
		self.imputers = [ SelfProgenyImputer(records, i, map_, 0.01)
							for i in range(self.num_progenies()) if i != ic ]
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 1
	
	def impute(self) -> None:
		self.parent_imputer.impute()
		for imputer in self.imputers:
			imputer.impute()
