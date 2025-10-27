from __future__ import annotations

# coding: utf-8

from collections import Counter
import sys

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from Genotype import Genotype
from TypeDeterminer import ParentComb


#################### VCFHomoHomoRecord ####################

class VCFHomoHomoRecord(VCFImpFamilyRecord):
	def __init__(self, pos: int, geno: list[int], index: int,
						parents_wrong_type: str, pair: ParentComb) -> None:
		super().__init__(pos, geno, index, parents_wrong_type, pair)
	
	def is_imputable(self) -> bool:
		return False
	
	def get_fill_type(self) -> FillType:
		if self.pair == ParentComb.P00x11 and (self.is_NA(0) or self.is_NA(1)):
			return FillType.IMPUTABLE
		else:
			return FillType.FILLED
	
	def gts(self) -> list[int]:
		if self.pair == ParentComb.P00x00:
			return [Genotype.PH_00] * self.num_samples()	# 0|0 ...
		elif self.pair == ParentComb.P11x11:
			return [Genotype.PH_11] * self.num_samples()	# 1|1 ...
		else:
			if ((self.is_mat_ref_homo() and not self.is_pat_ref_homo()) or
					(not self.is_mat_alt_homo() and self.is_pat_alt_homo())):
				# 0|0 1|1 0|1 0|1 ...
				return ([Genotype.PH_00, Genotype.PH_11] +
							[Genotype.PH_01] * self.num_progenies())
			if ((self.is_mat_alt_homo() and not self.is_pat_alt_homo()) or
					(not self.is_mat_ref_homo() and self.is_pat_ref_homo())):
				# 1|1 0|0 1|0 1|0 ...
				return ([Genotype.PH_11, Genotype.PH_00] +
							[Genotype.PH_10] * self.num_progenies())
			else:
				# ./. ./. 0/1 0/1 ...
				return ([Genotype.NA, Genotype.NA] +
							[Genotype.UN_01] * self.num_progenies())
	
	def impute(self) -> VCFHomoHomoRecord:
		return VCFHomoHomoRecord(self.pos, self.gts(), self.index,
											self.parents_wrong_type, self.pair)


#################### VCFHomoHomo ####################

class VCFHomoHomo(VCFFamilyBase):
	def __init__(self, samples: list[str],
					records: list[VCFHomoHomoRecord], vcf: VCFSmall) -> None:
		VCFFamilyBase.__init__(self, samples, vcf)
		self.records: list[VCFHomoHomoRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> list[VCFHomoHomo]:
		L: int = self.num_samples()
		records = [ record.impute() for record in self.records ]
		# HeteroHeteroがlistを返すのでそれに合わせる
		return [VCFHomoHomo(self.samples, records, self.vcf)]
