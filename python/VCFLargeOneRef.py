from __future__ import annotations

# coding: utf-8
# VCFLargeOneRef.py
# 片親がRefにあって、もう片親と後代がphasedとN/A

from collections import defaultdict, Counter
from math import log
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImputable import VCFImputable
from ParentRefImputer import *
from ProgenyRefImputer import *
from Map import *
from Genotype import Genotype


#################### VCFLargeOneRef ####################

class VCFLargeOneRef(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
								map_: Map, is_mat_ref: bool, ref_vcf: VCFGeno):
		VCFImputable.__init__(self, samples, ref_vcf.vcf)
		self.records: list[VCFFamilyRecord] = records
		self.is_mat_ref: bool = is_mat_ref
		self.parent_imputer = ParentRefImputer(records, map_,
												is_mat_ref, ref_vcf, 0.01)
		self.prog_imputer = ProgenyRefImputer(records, map_, 0.01)
	
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
	
	##### virtual methods for VCFImputable #####
	def impute(self) -> None:
		self.parent_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
