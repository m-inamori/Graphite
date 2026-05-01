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
	
	def copy(self) -> GenoRecord:
		return GenoRecord(self.pos, self.geno[:])
	
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
	
	@staticmethod
	def default_info(record: VCFRecord) -> str:
		w = record.v[9].split(':')
		w1 = []
		for u in w[1:]:
			t = u.split(',')
			u1 = ','.join(['.']*len(t))
			w1.append(u1)
		return ':'.join(w1)
	
	@staticmethod
	def extra_info(record: VCFRecord, c: int, def_info: str) -> str:
		if c == 0:
			return def_info
		else:
			return record.v[c][4:]
	
	# Write a VCF record to the output stream.
	#
	# This method outputs one VCF record, including per-sample fields.
	# For each sample, the genotype (GT) is always written based on the
	# internally stored integer representation.
	#
	# For additional per-sample FORMAT fields (i.e., fields other than GT):
	# - If the original extra information exists in the input record,
	#   it is preserved and written.
	# - Otherwise, default values are generated and used instead.
	#
	# The final record is written as a tab-separated line in standard VCF format.
	def write(self, record: VCFRecord, columns: list[int], out: TextIO) -> None:
		# Take the first 9 standard VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
		# from the given VCFRecord.
		v = record.v[:9]
		
		# Prepare the default per-sample auxiliary information (fields other than GT).
		# This will be used when the original extra information is not available.
		def_info = GenoRecord.default_info(record)
		
		# Iterate over all genotype entries stored in self.geno (one per sample).
		for i in range(len(self.geno)):
			# Get the integer-encoded genotype for the i-th sample.
			gt = self.geno[i]
			
			# Try to retrieve the original extra information (non-GT fields) for this sample.
			# If not available, fall back to the default information prepared above.
			ex_info = GenoRecord.extra_info(record, columns[i], def_info)
			
			# Convert the integer genotype to its string representation (e.g., "0/1").
			# If extra information exists, append it after ":" as in standard VCF FORMAT fields.
			if ex_info != '':
				s = Genotype.int_to_all_gt(gt) + ':' + ex_info
			else:
				s = Genotype.int_to_all_gt(gt)
			
			# Append the formatted sample field to the VCF record fields.
			v.append(s)
		
		# Write the complete VCF record as a tab-separated line to the output stream.
		write_tsv(v, out)
