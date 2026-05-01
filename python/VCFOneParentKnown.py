from __future__ import annotations

# coding: utf-8
# VCFOneParentKnown.py
# 片側だけ知られている


from collections import defaultdict, Counter
from math import log
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImputable import VCFImputable
from ParentImputer import *
from Map import *
from Genotype import Genotype


#################### VCFOneParentKnown ####################

class VCFOneParentKnown(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]],
					is_mat_known: bool, map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, vcf)
		self.records = records
		self.parent_imputer = ParentImputer(records, is_mat_known,
													ref_haps, map_, 0.01)
	
	##### virtual methods for VCFGenoBase #####
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_records(self) -> list[GenoRecord]:
		return [ r for r in self.records ]
	
	##### virtual methods for VCFFamilyBase #####
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	##### non-virtual methods #####
	def impute(self) -> None:
		self.parent_imputer.impute()
