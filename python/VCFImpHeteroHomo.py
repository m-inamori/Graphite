from __future__ import annotations

# coding: utf-8
# VCFImpHeteroHomo.py
# ヘテロ親がimputeされているヘテロ×ホモファミリー

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from VCFFillable import *
from group import Groups
from RecordSet import RecordSet
from VCFImpFamily import FillType
from VCFHeteroHomoOnePhased import *
from Map import *
import ClassifyRecord
import Imputer
from option import *


#################### VCFImpHeteroHomo ####################

class VCFImpHeteroHomo(VCFHeteroHomoOnePhased):
	def __init__(self, header: list[list[str]],
						records: list[VCFFillableRecord],
						is_mat_hetero: bool, map_: Map):
		VCFHeteroHomoOnePhased.__init__(self, header, records,
												is_mat_hetero, map_)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def imputed_index(self):
		return 0 if self.is_mat_hetero else 1
	
	def unimputed_index(self):
		return 1 if self.is_mat_hetero else 0
	
	def update_each(self, i: int, j: int, c: str) -> str:
		v = self.records[i].v
		k = int(c) * 2
		if self.is_mat_hetero:
			return v[9][k] + v[10][1:]
		else:
			return v[9][:2] + v[10][k]
	
	def update(self, i: int, seqs: list[str]):
		gt = self.records[i].get_gt(self.unimputed_index())
		self.records[i].set_GT(self.unimputed_index(), gt[0] + '|' + gt[2])
		for j in range(2, len(self.samples)):
			self.records[i].v[j+9] = self.update_each(i, j, seqs[j-2][i])
	
	# 0|1 1/1 0/0 -> '0'
	# 0|1 1/1 0/1 -> '0'
	# 1|0 1/1 0/1 -> '1'
	# 1|0 1/1 1/1 -> '0'
	def determine_haplotype(self, which_zero: int,
								homo_int_gt: int, prog_int_gt: int) -> str:
		if prog_int_gt == -1:
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
			gt = record.get_int_gt(i)
			if self.is_mat_hetero:
				which_zero = 0 if record.mat_gt()[0] == '0' else 1
				pat_int_gt = record.pat_int_gt()
				c = self.determine_haplotype(which_zero, pat_int_gt, gt)
			else:
				which_zero = 0 if record.pat_gt()[0] == '0' else 1
				mat_int_gt = record.mat_int_gt()
				c = self.determine_haplotype(which_zero, mat_int_gt, gt)
			cs.append(c)
		return ''.join(cs)
	
	def impute_sample_seq(self, j: int, cMs: list[float], min_c: float):
		seq = self.make_seq(j)
		if Imputer.is_all_same_without_N(seq):
			return Imputer.create_same_color_string(seq, '0')
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def determine_gts_from_unimputed_parent(self, j: int, hap: list[int]):
		j_hetero = 0 if self.is_mat_hetero else 1
		j_homo = 1 if self.is_mat_hetero else 0
		for i in range(len(self)):
			record = self.records[i]
			h = hap[i]
			gt_hetero = record.get_gt(j_hetero)[h*2]
			gt_homo = record.get_gt(j_homo)[0]
			if self.is_mat_hetero:
				record.set_GT(j, gt_hetero + '|' + gt_homo)
			else:
				record.set_GT(j, gt_homo + '|' + gt_hetero)
	
	def impute(self):
		if not self.records:
			return
		cMs = [ self.cM(record.pos()) for record in self.records ]
		imputed_seqs = [
				self.impute_sample_seq(i, cMs, 1.0)
								for i in range(2, self.num_samples()) ]
		for i in range(len(self)):
			self.update(i, imputed_seqs)
