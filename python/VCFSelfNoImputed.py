from __future__ import annotations

# coding: utf-8
# VCFSelfNoImputed.py
# 自殖でimputedなサンプルが一つもない

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from SelfImputer import SelfImputer
from Map import *


#################### VCFSelfNoImputed ####################

class VCFSelfNoImputed(VCFGenoBase, VCFMeasurable):
	def __init__(self, samples: list[str],
						records: list[GenoRecord],
						ref_haps: list[list[int]],
						map_: Map, vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records: list[GenoRecord] = records
		self.imputer = SelfImputer(records, ref_haps, map_, 0.01)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 1
	
	def impute(self) -> None:
		self.imputer.impute()
