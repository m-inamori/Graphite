from __future__ import annotations

# coding: utf-8
# VCFOneParentKnown.py
# 片側だけ知られている


from collections import defaultdict, Counter
from math import log
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from ParentImputer import *
from ProgenyImputer import *
from Map import *
from Genotype import Genotype


#################### VCFOneParentKnown ####################

class VCFOneParentKnown(VCFBase, VCFSmallBase):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_known: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		self.records = records
		self.parent1_imputer = ParentImputer(records, is_mat_known,
														ref_haps, map_, 0.01)
		self.parent2_imputer = ParentImputer(records, not is_mat_known,
														ref_haps, map_, 0.01)
		self.prog_imputer = ProgenyImputer(records, map_, 0.01)
	
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
	
	def impute_known_parent(self) -> None:
		self.parent1_imputer.impute()
	
	def impute(self) -> None:
		self.parent1_imputer.impute()
		self.parent2_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
