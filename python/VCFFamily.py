from __future__ import annotations

# coding: utf-8
# VCFFamily.py

from itertools import count, groupby
from abc import abstractmethod

from VCF import VCFRecord, VCFSmallBase, VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from pedigree import *
from Genotype import Genotype


#################### VCFFamilyRecord ####################

class VCFFamilyRecord(GenoRecord):
	def __init__(self, pos: int, geno: list[int]):
		super().__init__(pos, geno)
	
	def mat_gt(self) -> int:
		return self.geno[0]
	
	def pat_gt(self) -> int:
		return self.geno[1]
	
	def unphased_mat(self) -> int:
		return self.unphased(0)
	
	def unphased_pat(self) -> int:
		return self.unphased(1)
	
	def is_mat_ref_homo(self) -> bool:
		return self.is_ref_homo(0)
	
	def is_pat_ref_homo(self) -> bool:
		return self.is_ref_homo(1)
	
	def is_mat_alt_homo(self) -> bool:
		return self.is_alt_homo(0)
	
	def is_pat_alt_homo(self) -> bool:
		return self.is_alt_homo(1)
	
	def is_mat_hetero(self) -> bool:
		return self.is_hetero(0)
	
	def is_pat_hetero(self) -> bool:
		return self.is_hetero(1)
	
	def is_mat_NA(self) -> bool:
		return self.is_NA(0)
	
	def is_pat_NA(self) -> bool:
		return self.is_NA(1)
	
	def get_mat_allele(self, j: int) -> int:
		return self.get_allele(0, j)
	
	def get_pat_allele(self, j: int) -> int:
		return self.get_allele(1, j)
	
	def progeny_gts(self) -> list[int]:
		return self.geno[2:]
	
	def num_progenies(self) -> int:
		return len(self.geno) - 2


#################### VCFFamilyBase ####################

# Inherite this not VCFFamily for A family VCF
class VCFFamilyBase(VCFGenoBase):
	def __init__(self, samples: list[str], vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
	
	@abstractmethod
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		pass
	
	def mat(self) -> str:
		return self.get_samples()[0]
	
	def pat(self) -> str:
		return self.get_samples()[1]
	
	def parents(self) -> tuple[str, str]:
		return (self.mat(), self.pat())
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2


#################### VCFFamily ####################

class VCFFamily(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
														vcf: VCFSmall) -> None:
		VCFFamilyBase.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	@staticmethod
	def create(vcf: VCFSmall, samples: list[str]) -> VCFFamily:
		columns = vcf.extract_columns(samples)
		# [VCFFamilyRecord]
		records = [ VCFFamily.subset(r, columns) for r in vcf.records ]
		pos = [ int(r.v[1]) for r in vcf.records ]
		return VCFFamily(samples, records, vcf)
	
	# vcf1にあるsampleはvcf1から、
	# そうでないsampleはvcf2からGenotypeを取って新たなVCFを作る
	@staticmethod
	def create_by_two_vcfs(vcf1: VCFGenoBase, vcf2: VCFSmall,
										samples: list[str]) -> VCFFamily:
		# samplesは[mat, pat, prog1, ...]の前提
		columns1 = vcf1.extract_columns(samples)
		columns2 = vcf2.extract_columns(samples)
		new_records: list[VCFFamilyRecord] = []
		for i in range(len(vcf1)):
			record1: GenoRecord = vcf1.get_record(i)
			record2: VCFRecord = vcf2.get_record(i)
			geno: list[int] = []
			for i in range(len(samples)):
				if columns1[i] != -1:
					geno.append(record1.geno[columns1[i]])
				elif columns2[i] != -1:
					geno.append(Genotype.all_gt_to_int(record2.v[columns2[i]]))
				else:	# どちらのVCFにもそのサンプルは無い
					geno.append(Genotype.NA)
			new_record = VCFFamilyRecord(record1.pos, geno)
			new_records.append(new_record)
		return VCFFamily(samples, new_records, vcf2)
	
	@staticmethod
	def subset(record: VCFRecord, columns: list[int]) -> VCFFamilyRecord:
		v = record.v
		pos = int(v[1])
		geno = [ Genotype.all_gt_to_int(v[c]) if c != -1 else Genotype.NA
															for c in columns ]
		return VCFFamilyRecord(pos, geno)
	
	@staticmethod
	def convert(vcf: VCFFamilyBase) -> VCFFamily:
		records: list[VCFFamilyRecord] = [ vcf.get_family_record(i)
													for i in range(len(vcf)) ]
		return VCFFamily(vcf.get_samples(), records, vcf.get_ref_vcf())
