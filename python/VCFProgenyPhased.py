from __future__ import annotations

# coding: utf-8
# VCFProgenyPhased.py
# 子どもが一つでもphasingされた家系

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from Map import *
import Imputer


#################### VCFProgenyPhased ####################

class VCFProgenyPhased(VCFFamilyBase):
	def __init__(self, header: list[list[str]],
						records: list[VCFFamilyRecord], ppi: list[int]):
		VCFFamilyBase.__init__(self, header)
		self.records: list[VCFFamilyRecord] = records
		self.phased_progeny_indices: list[int] = ppi
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def determine_parent(self, is_mat):
		i = self.phased_progeny_indices[0]	# とりあえず先頭のサンプルを使う
		j = 0 if is_mat else 2
		c = 9 if is_mat else 10
		for k in range(len(self)):
			record = self.get_family_record(k)
			prog_GT = record.get_GT(i)
			int_gt = record.get_int_gt(c-9)
			remain = int_gt - int(prog_GT[j])
			if remain <= 0:
				record.v[c] = prog_GT[j] + '|0'
			else:
				record.v[c] = prog_GT[j] + '|1'
	
	def impute(self):
		if not self.records:
			return
		
		self.determine_parent(True)
		self.determine_parent(False)
	
	@staticmethod
	def impute_by_progeny(orig_vcf: VCFSmall, imputed_vcf: VCFSmall,
						samples: list[str], ppi: list[int]) -> VCFProgenyPhased:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		new_vcf = VCFProgenyPhased(vcf.header, vcf.records, ppi)
		new_vcf.impute()
		return new_vcf

__all__ = ['VCFProgenyPhased']
