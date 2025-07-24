from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCFFamily import *
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from TypeDeterminer import ParentComb


#################### VCFSelfHomoRecord ####################

class VCFSelfHomoRecord(VCFImpSelfRecord):
	def __init__(self, v: list[str], samples: list[str],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(v, samples, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> SelfFillType:
		return SelfFillType.FILLED
	
	def gts(self) -> list[str]:
		if self.pair == ParentComb.P00x00:
			return ['0|0'] * len(self.samples)
		else:
			return ['1|1'] * len(self.samples)
	
	def impute(self) -> VCFSelfHomoRecord:
		v = self.v[:9] + [ gt + s[3:] for gt, s in zip(self.gts(), self.v[9:]) ]
		return VCFSelfHomoRecord(v, self.samples, self.index,
											self.parents_wrong_type, self.pair)
