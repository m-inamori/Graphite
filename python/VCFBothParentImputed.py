from __future__ import annotations

# coding: utf-8
# VCFBothParentImputed.py
# 両親ともimputeされている

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from ProgenyImputer import *
from Genotype import Genotype


#################### VCFBothParentImputed ####################

class VCFBothParentImputed(VCFFamily):
	def __init__(self, samples: list[str],
					records: list[VCFFamilyRecord], map_: Map, vcf: VCFSmall):
		VCFFamily.__init__(self, samples, records, vcf)
		self.prog_imputer = ProgenyImputer(records, map_, 0.01)
		self.records: list[VCFFamilyRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_progenies(self) -> int:
		return len(self.get_samples()) - 2
	
	def impute(self, precision_ratio: float, num_threads: int) -> None:
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
