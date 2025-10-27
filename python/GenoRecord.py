from __future__ import annotations

# coding: utf-8
# GenoRecord.py
# 整数でGenotypeを持つ

from abc import ABC, abstractmethod
from typing import TextIO

from VCF import VCFRecord, VCFSmall
from Genotype import Genotype
from common import write_tsv


#################### GenoRecord ####################

class GenoRecord(ABC):
	def __init__(self, pos: int, geno: list[int]) -> None:
		self.pos: int = pos
		self.geno: list[int] = geno
	
	def num_samples(self) -> int:
		return len(self.geno)
	
	def unphased(self, i: int) -> int:
		return Genotype.unphased(self.geno[i])
	
	def is_ref_homo(self, i: int) -> bool:
		return Genotype.is_ref_homo(self.geno[i])
	
	def is_alt_homo(self, i: int) -> bool:
		return Genotype.is_alt_homo(self.geno[i])
	
	def is_hetero(self, i: int) -> bool:
		return Genotype.is_hetero(self.geno[i])
	
	def is_NA(self, i: int) -> bool:
		return Genotype.is_NA(self.geno[i])
	
	def is_00(self, i: int) -> bool:
		return Genotype.is_00(self.geno[i])
	
	def is_01(self, i: int) -> bool:
		return Genotype.is_01(self.geno[i])
	
	def is_11(self, i: int) -> bool:
		return Genotype.is_11(self.geno[i])
	
	def is_homo(self, i: int) -> bool:
		return Genotype.is_homo(self.geno[i])
	
	def is_phased(self, i: int) -> bool:
		return Genotype.is_phased(self.geno[i])
	
	# j(0 or 1)側のallele
	def get_allele(self, i: int, j: int) -> int:
		return Genotype.get_allele(self.geno[i], j)
	
	def unphased_gts(self) -> list[int]:
		return [ self.unphased(i) for i in range(self.num_samples()) ]
	
	def write(self, record: VCFRecord, out: TextIO) -> None:
		v = record.v[:9]
		for i in range(len(self.geno)):
			gt = self.geno[i]
			gt_orig = record.v[i+9]
			v.append(Genotype.int_to_all_gt(gt) + gt_orig[3:])
		write_tsv(v, out)
