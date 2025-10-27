from __future__ import annotations

# coding: utf-8
# VCFImpHeteroHomo.py
# ヘテロ親がimputeされているヘテロ×ホモファミリー

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpFamilyRecord import FillType
from VCFHeteroHomoOnePhased import *
from Map import *
from Genotype import Genotype
import Imputer
from option import *


#################### VCFImpHeteroHomo ####################

class VCFImpHeteroHomo(VCFHeteroHomoOnePhased):
	def __init__(self, samples: list[str],
						records: list[VCFFillableRecord],
						is_mat_hetero: bool, map_: Map, vcf: VCFSmall):
		VCFHeteroHomoOnePhased.__init__(self, samples, records,
												is_mat_hetero, map_, vcf)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def imputed_index(self) -> int:
		return 0 if self.is_mat_hetero else 1
	
	def unimputed_index(self) -> int:
		return 1 if self.is_mat_hetero else 0
	
	def update_each(self, i: int, j: int, c: str) -> int:
		record = self.records[i]
		k = int(c)
		if self.is_mat_hetero:
			a1 = record.get_allele(0, k)
			a2 = record.unphased(1) // 2
			return Genotype.from_alleles(a1, a2)
		else:
			a1 = record.unphased(0) // 2
			a2 = record.get_allele(1, k)
			return Genotype.from_alleles(a1, a2)
	
	def update(self, i: int, seqs: list[str]) -> None:
		record = self.records[i]
		a = record.unphased(self.unimputed_index()) // 2
		record.geno[self.unimputed_index()] = Genotype.from_alleles(a, a)
		for j in range(2, len(self.get_samples())):
			record.geno[j] = self.update_each(i, j, seqs[j-2][i])
	
	# 0|1 1/1 0/0 -> '0'
	# 0|1 1/1 0/1 -> '0'
	# 1|0 1/1 0/1 -> '1'
	# 1|0 1/1 1/1 -> '0'
	def determine_haplotype(self, which_zero: int,
								homo_int_gt: int, prog_int_gt: int) -> str:
		if prog_int_gt == Genotype.NA:
			return 'N'
		
		pat_int = homo_int_gt // 2
		for i in range(2):
			mat_int = (i + which_zero) & 1
			if mat_int + pat_int == prog_int_gt:
				return str(i)
			elif mat_int == 1 and mat_int + pat_int < prog_int_gt:
				return str(i)
			elif mat_int == 0 and mat_int + pat_int > prog_int_gt:
				return str(i)
		else:
			return 'N'
	
	def make_seq(self, i: int) -> str:
		cs: list[str] = []
		for record in self.records:
			gt = record.unphased(i)
			if self.is_mat_hetero:
				which_zero = record.get_mat_allele(0)
				pat_int_gt = record.unphased_pat()
				c = self.determine_haplotype(which_zero, pat_int_gt, gt)
			else:
				which_zero = record.get_pat_allele(0)
				mat_int_gt = record.unphased_mat()
				c = self.determine_haplotype(which_zero, mat_int_gt, gt)
			cs.append(c)
		return ''.join(cs)
	
	def impute_sample_seq(self, j: int, cMs: list[float], min_c: float) -> str:
		seq = self.make_seq(j)
		if Imputer.is_all_same_without_N(seq):
			return Imputer.create_same_color_string(seq, '0')
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def determine_gts_from_unimputed_parent(self, j: int,
												hap: list[int]) -> None:
		j_hetero = 0 if self.is_mat_hetero else 1
		j_homo = 1 if self.is_mat_hetero else 0
		for i in range(len(self)):
			record = self.records[i]
			h = hap[i]
			a_hetero = record.get_allele(j_hetero, h)
			a_homo = record.get_allele(j_homo, 0)
			if self.is_mat_hetero:
				record.geno[j] = Genotype.from_alleles(a_hetero, a_homo)
			else:
				record.geno[j] = Genotype.from_alleles(a_homo, a_hetero)
	
	def impute(self) -> None:
		if not self.records:
			return
		cMs = [ self.cM(record.pos) for record in self.records ]
		imputed_seqs = [
				self.impute_sample_seq(i, cMs, 1.0)
								for i in range(2, self.num_samples()) ]
		for i in range(len(self)):
			self.update(i, imputed_seqs)
