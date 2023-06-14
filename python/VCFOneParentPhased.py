from __future__ import annotations

# coding: utf-8
# VCFOneParentPhased.py
# 片方の親だけがphasingされた家系

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from VCFFillable import *
from Map import *
import Imputer
from option import *


#################### VCFHOneParentPhased ####################

class VCFOneParentPhased(VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]],
						records: list[VCFFamilyRecord], mat_p: bool, map_: Map):
		VCFFamilyBase.__init__(self, header)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.is_mat_phased: bool = mat_p
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def determine_which_comes_from(self, record: VCFFamilyRecord,
														i: int) -> str:
		if record.is_NA(0) or record.is_NA(1) or record.is_NA(i):
			return 'N'
		
		mat_GT = record.get_GT(0)
		pat_GT = record.get_GT(1)
		prog_gt = record.get_int_gt(i)
		mat1 = int(mat_GT[0])
		mat2 = int(mat_GT[2])
		pat1 = int(pat_GT[0])
		pat2 = int(pat_GT[2])
		if self.is_mat_phased:
			if not (mat1 != mat2 and pat1 == pat2):
				return 'N'
			
			diff = prog_gt - pat1
			if diff == -1 or diff == 2:
				return 'N'
			else:
				return '0' if diff == mat1 else '1'
		else:	# pat is phased
			if not (mat1 == mat2 and pat1 != pat2):
				return 'N'
			
			diff = prog_gt - mat1
			if diff == -1 or diff == 2:
				return 'N'
			else:
				return '0' if diff == pat1 else '1'
	
	def make_seq(self, i: int) -> str:
		cs = []
		for record in self.records:
			cs.append(self.determine_which_comes_from(record, i))
		
		return ''.join(cs)
	
	def impute_sample_seq(self, i: int, cMs: list[float], min_c: float):
		seq = self.make_seq(i)
		if all(c == 'N' for c in seq):
			# 全部同じなら意味がない
			return '0' * len(seq)
		elif is_all_same(seq):
			return seq
		
		hidden_states = ['0', '1']
		states = [ s for s in ['0', '1', 'N'] if s in seq ]
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def update_each(self, i: int, k: int, c: str):
		record = self.records[k]
		v = record.v
		j = int(c) * 2
		if self.is_mat_phased:
			mat_gt = v[9][j]
			remain = record.get_int_gt(i) - int(mat_gt)
			if remain <= 0:
				record.set_GT(i, mat_gt + '|0')
			else:
				record.set_GT(i, mat_gt + '|1')
		else:
			pat_gt = v[10][j]
			remain = record.get_int_gt(i) - int(pat_gt)
			if remain <= 0:
				record.set_GT(i, '0|' + pat_gt)
			else:
				record.set_GT(i, '1|' + pat_gt)
		
		if '0' in self.samples:
			return
		
		if self.is_mat_phased:
			mat_gt = v[9][j]
			remain = record.get_int_gt(i) - int(mat_gt)
			if remain <= 0:
				record.set_GT(1, '0|0' if record.get_int_gt(1) == 0 else '0|1')
			else:
				record.set_GT(1, '1|1' if record.get_int_gt(1) == 2 else '1|0')
		else:
			pat_gt = v[10][j]
			remain = record.get_int_gt(i) - int(pat_gt)
			if remain <= 0:
				record.set_GT(0, '0|0' if record.get_int_gt(0) == 0 else '0|1')
			else:
				record.set_GT(0, '1|1' if record.get_int_gt(0) == 2 else '1|0')
	
	def update(self, i: int, seqs: str):
		for k in range(len(self)):
			self.update_each(i, k, seqs[k])
	
	def impute(self):
		if not self.records:
			return
		
		cMs = [ self.cM(record.pos()) for record in self.records ]
		for i in range(2, self.num_samples()):
			imputed_seq = self.impute_sample_seq(i, cMs, 1.0)
			self.update(i, imputed_seq)
	
	@staticmethod
	def impute_by_parent(orig_vcf: VCFSmall, imputed_vcf: VCFSmall,
										samples: list[str], is_mat_phased: bool,
										gmap: Map) -> VCFOneParentPhased:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		new_vcf = VCFOneParentPhased(vcf.header, vcf.records,
														is_mat_phased, gmap)
		new_vcf.impute()
		return new_vcf
	
	# samplesに対応するRecordを作る
	@staticmethod
	def merge_records(vcfs: list[VCFFamily], i: int,
									samples: list[str]) -> VCFRecord:
		dic_pos: dict[str, tuple[VCFFamily, int]] = { }
		for vcf in vcfs:
			for k, sample in enumerate(vcf.samples):
				dic_pos[sample] = (vcf, k + 9)
		
		v = vcfs[0].records[i].v[:9]
		for sample in samples:
			vcf, k = dic_pos[sample]
			v.append(vcf.records[i].v[k])
		return VCFRecord(v, samples)

__all__ = ['VCFOneParentPhased']
