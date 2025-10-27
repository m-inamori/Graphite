from __future__ import annotations

# coding: utf-8
# imputeしないバージョン
# selfのときのみimputeするバージョンを使う

from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb


#################### VCFHeteroHeteroLiteRecord ####################

class VCFHeteroHeteroLiteRecord(VCFImpFamilyRecord):
	def __init__(self, pos: int, geno: list[int],
							index: int, parents_wrong_type: str) -> None:
		super().__init__(pos, geno, index,
								parents_wrong_type, ParentComb.P01x01)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		return FillType.IMPUTABLE
