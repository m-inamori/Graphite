from __future__ import annotations

# coding: utf-8
# VCFJunkRecord.py

from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from TypeDeterminer import ParentComb


#################### VCFJunkRecord ####################

class VCFSelfJunkRecord(VCFImpSelfRecord):
	def __init__(self, pos: int, geno: list[int],
								index: int, parents_wrong_type: str):
		super().__init__(pos, geno, index, parents_wrong_type, ParentComb.PNA)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> SelfFillType:
		return SelfFillType.UNABLE
