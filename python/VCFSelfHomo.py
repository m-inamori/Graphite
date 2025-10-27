from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from Genotype import Genotype
from TypeDeterminer import ParentComb


#################### VCFSelfHomoRecord ####################

class VCFSelfHomoRecord(VCFImpFamilyRecord):
	def __init__(self, pos: int, geno: list[int],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(pos, geno, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		return FillType.FILLED
	
	def gts(self) -> list[int]:
		if self.pair == ParentComb.P00x00:
			return [Genotype.PH_00] * len(self.geno)
		else:
			return [Genotype.PH_11] * len(self.geno)
	
	def impute(self) -> VCFSelfHomoRecord:
		return VCFSelfHomoRecord(self.pos, self.gts(), self.index,
											self.parents_wrong_type, self.pair)


#################### VCFSelfHomo ####################

class VCFSelfHomo(VCFGenoBase):
	def __init__(self, samples: list[str],
						records: list[VCFSelfHomoRecord], vcf: VCFSmall):
		VCFGenoBase.__init__(self, samples, vcf)
		self.records: list[VCFSelfHomoRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> list[VCFSelfHomo]:
		L: int = self.num_samples()
		records = [ record.impute() for record in self.records ]
		# HeteroHeteroがlistを返すのでそれに合わせる
		return [VCFSelfHomo(self.samples, records, self.vcf)]
