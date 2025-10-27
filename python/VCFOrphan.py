from __future__ import annotations

# coding: utf-8
# VCFOrphan.py
# 孤立したサンプルを集めたVCF

from typing import List, Tuple

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from OrphanImputer import *
from Genotype import Genotype


#################### VCFOrphan ####################

DP = List[Tuple[float, int]]	# (log of probability, prev h)
MIN_PROB = -1e300

class VCFOrphan(VCFGenoBase):
	def __init__(self, samples: list[str], records: list[GenoRecord],
						ref_haps: list[list[int]], map_: Map, vcf: VCFSmall):
		VCFGenoBase.__init__(self, samples, vcf)
		self.imputer = OrphanImputer(records, ref_haps, map_, 0.01)
		self.records: list[GenoRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def impute(self) -> None:
		for i in range(self.num_samples()):
			self.imputer.impute(i)
