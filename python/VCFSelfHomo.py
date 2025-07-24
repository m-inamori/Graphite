from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCFFamily import *
from VCFImpFamily import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb


#################### VCFSelfHomoRecord ####################

class VCFSelfHomoRecord(VCFImpFamilyRecord):
	def __init__(self, v: list[str], samples: list[str],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(v, samples, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		return FillType.FILLED
	
	def gts(self) -> list[str]:
		if self.pair == ParentComb.P00x00:
			return ['0|0'] * len(self.samples)
		else:
			return ['1|1'] * len(self.samples)
	
	def impute(self) -> VCFSelfHomoRecord:
		v = self.v[:9] + [ gt + s[3:] for gt, s in zip(self.gts(), self.v[9:]) ]
		return VCFSelfHomoRecord(v, self.samples, self.index,
											self.parents_wrong_type, self.pair)


#################### VCFSelfHomo ####################

class VCFSelfHomo(VCFBase, VCFSmallBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFSelfHomoRecord]):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		self.records = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> list[VCFSelfHomo]:
		L: int = self.num_samples()
		records = [ record.impute() for record in self.records ]
		# HeteroHeteroがlistを返すのでそれに合わせる
		return [VCFSelfHomo(self.header, records)]
