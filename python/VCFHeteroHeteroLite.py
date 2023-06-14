from __future__ import annotations

# coding: utf-8
# imputeしないバージョン
# selfのときのみimputeするバージョンを使う

from VCFFamily import *
from VCFImpFamily import VCFImpFamilyRecord


#################### VCFHeteroHeteroLiteRecord ####################

class VCFHeteroHeteroLiteRecord(VCFImpFamilyRecord):
	def __init__(self, v: list[str], samples: list[str],
								index: int, parents_wrong_type: str):
		super().__init__(v, samples, index, parents_wrong_type, 2)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> str:
		return 'IMPUTABLE'


#################### VCFHeteroHeteroLite ####################

class VCFHeteroHeteroLite(VCFFamilyBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFHeteroHeteroLiteRecord]):
		self.records = records
		VCFFamilyBase.__init__(self, header)
