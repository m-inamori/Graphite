from __future__ import annotations

# coding: utf-8
# VCFOrphan.py
# 孤立したサンプルを集めたVCF

from typing import List, Tuple

from VCFFamily import *
from OrphanImputer import *
from Genotype import Genotype


#################### VCFOrphan ####################

DP = List[Tuple[float, int]]	# (log of probability, prev h)
MIN_PROB = -1e300

class VCFOrphan(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFRecord],
								ref_haps: list[list[int]], map_: Map) -> None:
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.imputer = OrphanImputer(records, ref_haps, map_, 0.01)
		self.records: list[VCFRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def impute(self) -> None:
		for i in range(self.num_samples()):
			self.imputer.impute(i)
