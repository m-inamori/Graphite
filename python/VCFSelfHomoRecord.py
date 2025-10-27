from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCFFamily import *
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from Genotype import Genotype
from TypeDeterminer import ParentComb


#################### VCFSelfHomoRecord ####################

class VCFSelfHomoRecord(VCFImpSelfRecord):
	def __init__(self, pos: int, geno: list[int],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(pos, geno, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> SelfFillType:
		return SelfFillType.FILLED
	
	def gts(self) -> list[int]:
		if self.pair == ParentComb.P00x00:
			return [Genotype.PH_00] * self.num_samples()
		else:
			return [Genotype.PH_11] * self.num_samples()
	
	def impute(self) -> VCFSelfHomoRecord:
		return VCFSelfHomoRecord(self.pos, self.gts(), self.index,
											self.parents_wrong_type, self.pair)
