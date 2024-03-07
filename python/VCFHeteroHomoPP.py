from __future__ import annotations

# coding: utf-8
# VCFHeteroHomoPP.py
# 両親がphasingされたヘテロ×ホモファミリー

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from VCFFillable import *
from RecordSet import RecordSet
from VCFImpFamily import FillType
from Map import *
from TypeDeterminer import ParentComb, TypeDeterminer
import ClassifyRecord
import Imputer
from option import *


#################### VCFHeteroHomoPP ####################

class VCFHeteroHomoPP(VCFBase, VCFSmallBase, VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]],
						records: list[VCFFillableRecord], map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def make_seq(self, i):
		cs = []
		for record in self.records:
			mat_GT = record.get_GT(0)
			pat_GT = record.get_GT(1)
			prog_gt = record.get_int_gt(i+2)
			if prog_gt == -1:
				cs.append('N')
				continue
			
			mat1 = int(mat_GT[0])
			mat2 = int(mat_GT[2])
			pat1 = int(pat_GT[0])
			pat2 = int(pat_GT[2])
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
	
	def impute_sample_seq(self, i: int, cMs: list[float], min_c: float):
		seq = self.make_seq(i)
		if VCFHeteroHomoPP.is_all_same_without_N(seq):
			return VCFHeteroHomoPP.create_same_color_string(seq)
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def update_each(self, i: int, j: int, c: str) -> str:
		v = self.records[i].v
		k = int(c) * 2
		if v[9][0] != v[9][2]:		# mat hetero
			return v[9][k] + v[10][1:]
		else:
			return v[9][:2] + v[10][k]
	
	def update(self, i: int, seqs: list[str]):
		for j in range(2, len(self.samples)):
			self.records[i].v[j+9] = self.update_each(i, j, seqs[j-2][i])
	
	def impute(self):
		if not self.records:
			return
		cMs = [ self.cM(record.pos()) for record in self.records ]
		imputed_seqs = [
				self.impute_sample_seq(i, cMs, 1.0)
								for i in range(self.num_samples() - 2) ]
		for i in range(len(self)):
			self.update(i, imputed_seqs)
	
	def fill(self):
		# FillTypeでrecordを分ける
		groups: list[tuple[FillType, list[VCFFillableRecord]]] = [ (g, list(v))
					for g, v in groupby(self.records, key=lambda r: r.type) ]
		for i, (key, subrecords) in enumerate(groups):
			if key != FillType.MAT and key != FillType.PAT:
				self.__phase(groups, i, True)
	
	# このあたりはあとでVCFFillableと共通化する
	def __phase(self, groups: list[tuple[FillType, list[VCFFillableRecord]]],
											i: int, necessary_parents_phasing):
		prev_mat_record = self.find_prev_record(groups, i, FillType.MAT)
		next_mat_record = self.find_next_record(groups, i, FillType.MAT)
		prev_pat_record = self.find_prev_record(groups, i, FillType.PAT)
		next_pat_record = self.find_next_record(groups, i, FillType.PAT)
		key, records = groups[i]
		for record in records:
			record_set = RecordSet(record, prev_mat_record,
							next_mat_record, prev_pat_record, next_pat_record)
			self.__impute_core(record_set)
	
	def __impute_core(self, record_set: RecordSet):
		record = record_set.record
		if record is None:
			return
		for i in range(2, len(record.samples)):
			_, mat_gt1, mat_gt2, pat_gt1, pat_gt2 = record_set.gts(i)
			prev_mat_from = record_set.from_which_chrom_prev_mat(mat_gt1)
			next_mat_from = record_set.from_which_chrom_next_mat(mat_gt2)
			prev_pat_from = record_set.from_which_chrom_prev_pat(pat_gt1)
			next_pat_from = record_set.from_which_chrom_next_pat(pat_gt2)
			# とりあえず、両側Noneはないと仮定する
			if prev_mat_from == 0:
				mat_from = next_mat_from
			elif next_mat_from == 0:
				mat_from = prev_mat_from
			elif prev_mat_from == next_mat_from:
				mat_from = prev_mat_from
			elif record_set.is_prev_nearer(True):
				mat_from = prev_mat_from
			else:
				mat_from = next_mat_from
			if prev_pat_from == 0:
				pat_from = next_pat_from
			elif next_pat_from == 0:
				pat_from = prev_pat_from
			elif prev_pat_from == next_pat_from:
				pat_from = prev_pat_from
			elif record_set.is_prev_nearer(True):
				pat_from = prev_pat_from
			else:
				pat_from = next_pat_from
			v = record.v
			v[i+9] = v[9][mat_from*2-2] + '|' + v[10][pat_from*2-2]
	
	def find_prev_record(self,
						groups: list[tuple[FillType, list[VCFFillableRecord]]],
						i: int, g: FillType) -> Optional[VCFFillableRecord]:
		key, records = groups[i]
		chr = records[0].v[0]
		for j in range(i-1, -1, -1):
			key, records = groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[-1]
		else:
			return None
	
	def find_next_record(self,
						groups: list[tuple[FillType, list[VCFFillableRecord]]],
						i: int, g: FillType) -> Optional[VCFFillableRecord]:
		key, records = groups[i]
		chr = records[0].v[0]
		for j in range(i+1, len(groups)):
			key, records = groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[0]
		else:
			return None
	
	@staticmethod
	def classify_record(record: VCFRecord) -> tuple[ParentComb, FillType]:
		i = 0 if record.v[9][0] != record.v[9][2] else 1
		j = 0 if record.v[10][0] != record.v[10][2] else 1
		if (i, j) == (0, 0):		# 0/1 x 0/1
			return (ParentComb.P01x01, FillType.IMPUTABLE)
		elif (i, j) == (1, 0):		# 0/0 x 0/1 or 1/1 x 0/1
			if record.v[9][0] == '0':
				return (ParentComb.P00x01, FillType.PAT)
			else:
				return (ParentComb.P01x11, FillType.PAT)
		elif (i, j) == (0, 1):		# 0/1 x 0/0 or 0/1 x 1/1
			if record.v[10][0] == '0':
				return (ParentComb.P00x01, FillType.MAT)
			else:
				return (ParentComb.P01x11, FillType.MAT)
		else:
			if record.v[9][0] == '0' and record.v[10][0] == '0':	# 0/0 x 0/0
				return (ParentComb.P00x00, FillType.FILLED)
			elif record.v[9][0] == '1' and record.v[10][0] == '1':	# 1/1 x 1/1
				return (ParentComb.P11x11, FillType.FILLED)
			else:													# 0/0 x 1/1
				return (ParentComb.P00x11, FillType.FILLED)
	
	@staticmethod
	def classify_records(records: list[VCFFamilyRecord]
									) -> list[list[VCFFillableRecord]]:
		# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
		rss: list[list[VCFFillableRecord]] = [ [] for _ in range(4) ]
		for index, record in enumerate(records):
			pair, type = VCFHeteroHomoPP.classify_record(record)
			rss[type.value].append(VCFFillableRecord(record.v, record.samples,
															index, type, pair))
		return rss
	
	@staticmethod
	def merge_vcf(mat_vcf: VCFHeteroHomoPP, pat_vcf: VCFHeteroHomoPP,
					homohomo_records: list[VCFFillableRecord],
					heterohetero_records: list[VCFFillableRecord]
														) -> VCFFillable:
		for record in homohomo_records:
			GT = record.v[9][0] + '|' + record.v[10][0]
			for j in range(2, len(record.samples)):
				record.v[j+9] = GT
		
		records = sorted(mat_vcf.records + pat_vcf.records +
									homohomo_records + heterohetero_records,
														key=lambda r: r.pos())
		return VCFFillable(mat_vcf.header, records)
	
	@staticmethod
	def impute_by_parents(orig_vcf: VCFSmall, imputed_vcf: VCFSmallBase,
								samples: list[str], gmap: Map) -> VCFFillable:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
		rss = VCFHeteroHomoPP.classify_records(vcf.records)
		mat_vcf = VCFHeteroHomoPP(vcf.header, rss[FillType.MAT.value], gmap)
		pat_vcf = VCFHeteroHomoPP(vcf.header, rss[FillType.PAT.value], gmap)
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
		v = record1.v + record2.v[9:]
		record = VCFFamilyRecord(v, samples)
		pair, type = VCFHeteroHomoPP.classify_record(record)
		return VCFFillableRecord(v, samples, i, type, pair)
	
	@staticmethod
	def fill_NA(record1: VCFRecord,
						samples: list[str], i: int) -> VCFFillableRecord:
		NA_len = len(samples) - len(record1.samples)
		v = record1.v + ['./.'] * NA_len
		return VCFFillableRecord(v, samples, i, FillType.UNABLE, ParentComb.PNA)
	
	@staticmethod
	def merge(vcf_parents: VCFSmall, vcf_progenies: VCFSmall,
				samples: list[str], m: Map, option: Option) -> VCFHeteroHomoPP:
		td = ClassifyRecord.get_typedeterminer(len(samples)-2, option.ratio)
		header = vcf_parents.trim_header(samples)
		records: list[VCFFillableRecord] = []
		j = 0
		for i, record1 in enumerate(vcf_parents.records):
			if j == len(vcf_progenies):
				record = VCFHeteroHomoPP.fill_NA(record1, samples, i)
			else:
				record2 = vcf_progenies.records[j]
				if record1.pos() == record2.pos():
					record = VCFHeteroHomoPP.merge_record(record1, record2,
																samples, i, td)
					j += 1
				else:
					record = VCFHeteroHomoPP.fill_NA(record1, samples, i);
			records.append(record)
		return VCFHeteroHomoPP(header, records, m)

__all__ = ['VCFHeteroHomoPP']
