from __future__ import annotations

# coding: utf-8
# VCFNoParentImputed.py
# 両親ともimputeされていない
# 母親を補完する

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple

from VCFFamily import *
from ParentImputer import *
from ProgenyImputer import *
from Genotype import Genotype


#################### VCFNoParentImputed ####################

DP = List[Tuple[float, int]]	# (log of probability, prev h)
MIN_PROB = -1e300

class VCFNoParentImputed(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
										ref_haps: list[list[int]], map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.mat_imputer = ParentImputer(records, True, ref_haps, map_, 0.01)
		self.pat_imputer = ParentImputer(records, False, ref_haps, map_, 0.01)
		self.prog_imputer = ProgenyImputer(records, map_, 0.01)
		self.records: list[VCFFamilyRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2
	
	def impute(self):
		self.mat_imputer.impute()
		self.pat_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
