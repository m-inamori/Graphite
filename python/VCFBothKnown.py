from __future__ import annotations

# coding: utf-8
# VCFBothKnown.py
# 両親ともimputeされていない
# 母親を補完する

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from ParentImputer import *
from ProgenyImputer import *
from Genotype import Genotype


#################### VCFBothKnown ####################

DP = List[Tuple[float, int]]	# (log of probability, prev h)
MIN_PROB = -1e300

class VCFBothKnown(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
						ref_haps1: list[list[int]], ref_haps2: list[list[int]],
						map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		self.mat_imputer = ParentImputer(records, True, ref_haps1, map_, 0.01)
		self.pat_imputer = ParentImputer(records, False, ref_haps2, map_, 0.01)
		self.prog_imputer = ProgenyImputer(records, map_, 0.01)
		self.records: list[VCFFamilyRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> None:
		self.mat_imputer.impute()
		self.pat_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
