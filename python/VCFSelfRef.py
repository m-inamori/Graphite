from __future__ import annotations

# coding: utf-8
# VCFSelfRef.py
# 自殖で親がRefにあって後代がphasedとN/A

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFSelfImputable import *
from ProgenySelfRefImputer import *
from Map import *


#################### VCFSelfRef ####################

class VCFSelfRef(VCFSelfImputable):
	def __init__(self, samples: list[str], records: list[GenoRecord],
													map_: Map, vcf: VCFSmall):
		VCFSelfImputable.__init__(self, samples, vcf)
		self.records: list[GenoRecord] = records
		self.imputer = ProgenySelfRefImputer(records, map_, 0.01)
	
	##### virtual methods for VCFGenoBase #####
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_records(self) -> list[GenoRecord]:
		return self.records
	
	##### virtual methods for VCFSelfImputable #####
	def impute(self) -> None:
		for i in range(self.num_progenies()):
			self.imputer.impute(i)
