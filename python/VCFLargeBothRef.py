from __future__ import annotations

# coding: utf-8
# VCFLargeBothRef.py
# 両親がRefにあって、後代がphasedとN/A

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImputable import *
from ProgenyRefImputer import *
from Map import *


#################### VCFLargeBothRef ####################

class VCFLargeBothRef(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
													map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.imputer = ProgenyRefImputer(records, map_, 0.01)
	
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
		for i in range(self.num_progenies()):
			self.imputer.impute(i)
