from __future__ import annotations

# coding: utf-8
# VCFHeteroHomoPP.py
# 両親がphasingされたヘテロ×ホモファミリー

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFRecord, VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFFillableRecord import *
from VCFFillable import *
from group import Groups
from RecordSet import RecordSet
from VCFImpFamilyRecord import FillType
from Map import *
from TypeDeterminer import ParentComb, TypeDeterminer
import ClassifyRecord
import Imputer
from Genotype import Genotype
from option import *


#################### VCFHeteroHomoPP ####################

class VCFHeteroHomoPP(VCFFamilyBase, VCFMeasurable):
	def __init__(self, samples: list[str],
						records: list[VCFFillableRecord],
						map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def make_seq(self, i: int) -> str:
		cs: list[str] = []
		for record in self.records:
			if record.is_NA(i+2):
				cs.append('N')
				continue
			
			mat1 = record.get_allele(0, 0)
			mat2 = record.get_allele(0, 1)
			pat1 = record.get_allele(1, 0)
			pat2 = record.get_allele(1, 1)
			prog_gt = record.unphased(i+2)
			if pat1 == pat2:	# mat hetero
				diff = prog_gt - pat1
				if diff == -1 or diff == 2:
					cs.append('N')
				else:
					cs.append('0' if diff == mat1 else '1')
			else:				# pat hetero
				diff = prog_gt - mat1
				if diff == -1 or diff == 2:
					cs.append('N')
				else:
					cs.append('0' if diff == pat1 else '1')
		
		return ''.join(cs)
	
	@staticmethod
	def is_all_same_without_N(seq: str) -> bool:
		c: str = '.'
		for c1 in seq:
			if c1 != 'N':
				if c == '.':
					c = c1
				elif c1 != c:
					return False
		else:
			return True
	
	@staticmethod
	def create_same_color_string(seq: str) -> str:
		c = '0'		# dummy
		for c1 in seq:
			if c1 != 'N':
				c = c1
				break
		return c * len(seq)
	
	def impute_sample_seq(self, i: int, cMs: list[float], min_c: float) -> str:
		seq = self.make_seq(i)
		if VCFHeteroHomoPP.is_all_same_without_N(seq):
			return VCFHeteroHomoPP.create_same_color_string(seq)
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def update_each(self, i: int, c: str) -> int:
		record = self.records[i]
		k = int(c)
		if record.is_mat_hetero():
			a1 = record.get_allele(0, k)
			a2 = record.get_allele(1, 1)
			return a1 | (a2 << 1) | 4
		else:
			a1 = record.get_allele(0, 0)
			a2 = record.get_allele(1, k)
			return a1 | (a2 << 1) | 4
	
	def update(self, i: int, seqs: list[str]) -> None:
		for j in range(2, len(self.samples)):
			self.records[i].geno[j] = self.update_each(i, seqs[j-2][i])
	
	def impute(self) -> None:
		if not self.records:
			return
		cMs = [ self.cM(record.pos) for record in self.records ]
		imputed_seqs = [
				self.impute_sample_seq(i, cMs, 1.0)
								for i in range(self.num_samples() - 2) ]
		for i in range(len(self)):
			self.update(i, imputed_seqs)
	
	def fill(self) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets():
			self.__impute_core(record_set)
	
	def __impute_core(self, record_set: RecordSet) -> None:
		record = record_set.record
		if record is None:
			return
		for i in range(2, len(record.geno)):
			mat_from = record_set.determine_mat_from(i)
			pat_from = record_set.determine_pat_from(i)
			a1 = record.get_allele(0, mat_from - 1)
			a2 = record.get_allele(1, pat_from - 1)
			record.geno[i] = a1 | (a2 << 1) | 4
	
	@staticmethod
	def classify_record(record: VCFFamilyRecord) -> tuple[ParentComb, FillType]:
		i = 0 if record.is_mat_hetero() else 1
		j = 0 if record.is_pat_hetero() else 1
		if (i, j) == (0, 0):		# 0/1 x 0/1
			return (ParentComb.P01x01, FillType.IMPUTABLE)
		elif (i, j) == (1, 0):		# 0/0 x 0/1 or 1/1 x 0/1
			if record.is_00(0):
				return (ParentComb.P00x01, FillType.PAT)
			else:
				return (ParentComb.P01x11, FillType.PAT)
		elif (i, j) == (0, 1):		# 0/1 x 0/0 or 0/1 x 1/1
			if record.is_00(1):
				return (ParentComb.P00x01, FillType.MAT)
			else:
				return (ParentComb.P01x11, FillType.MAT)
		else:
			if record.is_00(0) and record.is_00(1):
				return (ParentComb.P00x00, FillType.FILLED)
			elif record.is_11(0) and record.is_11(1):
				return (ParentComb.P11x11, FillType.FILLED)
			else:													# 0/0 x 1/1
				return (ParentComb.P00x11, FillType.FILLED)
	
	@staticmethod
	def classify_records(samples: list[str],
							records: Sequence[VCFFamilyRecord],
							ref_vcf: VCFSmall) -> list[list[VCFFillableRecord]]:
		cols = ref_vcf.extract_columns(samples)
		# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
		rss: list[list[VCFFillableRecord]] = [ [] for _ in range(4) ]
		for index, record in enumerate(records):
			pair, type = VCFHeteroHomoPP.classify_record(record)
			ref_record = ref_vcf.records[index]
			probs = ref_record.parse_PL(record.geno, cols)
			rss[type.value].append(VCFFillableRecord(record.pos, record.geno,
													index, type, pair, probs))
		return rss
	
	@staticmethod
	def merge_vcf(mat_vcf: VCFHeteroHomoPP, pat_vcf: VCFHeteroHomoPP,
					homohomo_records: list[VCFFillableRecord],
					heterohetero_records: list[VCFFillableRecord]
														) -> VCFFillable:
		for record in homohomo_records:
			a1 = record.get_allele(0, 0)
			a2 = record.get_allele(1, 0)
			gt = a1 | (a2 << 1) | 4
			for j in range(2, len(record.geno)):
				record.geno[j] = gt
		
		records = sorted(mat_vcf.records + pat_vcf.records +
									homohomo_records + heterohetero_records,
														key=lambda r: r.pos)
		return VCFFillable(mat_vcf.samples, records, mat_vcf.vcf)
	
	@staticmethod
	def impute_by_parents(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
								samples: list[str], gmap: Map) -> VCFFillable:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
		rss = VCFHeteroHomoPP.classify_records(samples, vcf.records, orig_vcf)
		mat_vcf = VCFHeteroHomoPP(samples, rss[FillType.MAT.value],
															gmap, orig_vcf)
		pat_vcf = VCFHeteroHomoPP(samples, rss[FillType.PAT.value],
															gmap, orig_vcf)
		mat_vcf.impute()
		pat_vcf.impute()
		new_vcf = VCFHeteroHomoPP.merge_vcf(mat_vcf, pat_vcf,
											rss[FillType.FILLED.value],
											rss[FillType.IMPUTABLE.value])
		new_vcf.phase_hetero_hetero()
		return new_vcf
	
	@staticmethod
	def merge_record(record1: VCFRecord, record2: VCFRecord,
						samples: list[str],
						i: int, td: TypeDeterminer) -> VCFFillableRecord:
		# tdがclassifyに使われていないのがおかしい
		pos = int(record1.v[1])
		geno = [ Genotype.all_gt_to_int(gt)
							for gt in record1.v[9:] + record2.v[9:] ]
		record = VCFFamilyRecord(pos, geno)
		pair, type = VCFHeteroHomoPP.classify_record(record)
		probs = [ VCFRecord.decide_PL_by_genotype(gt) for gt in geno ]
		return VCFFillableRecord(pos, geno, i, type, pair, probs)
	
	@staticmethod
	def fill_NA(record1: VCFRecord, samples: list[str],
												i: int) -> VCFFillableRecord:
		pos = int(record1.v[1])
		NA_len = len(samples) - len(record1.samples)
		geno = ([ Genotype.all_gt_to_int(gt) for gt in record1.v[9:] ] +
														[Genotype.NA] * NA_len)
		probs = [ VCFRecord.decide_PL_by_genotype(gt) for gt in geno ]
		return VCFFillableRecord(pos, geno, i, FillType.UNABLE,
												ParentComb.PNA, probs)
	
	@staticmethod
	def merge(vcf_parents: VCFSmall, vcf_progenies: VCFSmall,
						orig_vcf: VCFSmall, samples: list[str],
						m: Map, option: Option) -> VCFHeteroHomoPP:
		# 後代で無いポジションはN/Aで埋める
		td = ClassifyRecord.get_typedeterminer(len(samples)-2, option.ratio)
		records: list[VCFFillableRecord] = []
		j = 0
		for i, record1 in enumerate(vcf_parents.records):
			ref_record = orig_vcf.records[i]
			if j == len(vcf_progenies):
				record = VCFHeteroHomoPP.fill_NA(record1, samples, i)
			else:
				record2 = vcf_progenies.records[j]
				if record1.pos == record2.pos:
					record = VCFHeteroHomoPP.merge_record(record1, record2,
																samples, i, td)
					j += 1
				else:
					record = VCFHeteroHomoPP.fill_NA(record1, samples, i);
			records.append(record)
		return VCFHeteroHomoPP(samples, records, m, orig_vcf)

__all__ = ['VCFHeteroHomoPP']
