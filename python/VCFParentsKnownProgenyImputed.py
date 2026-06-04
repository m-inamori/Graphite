from __future__ import annotations

# coding: utf-8
# VCFParentsKnownProgenyImputed.py
# 両親がknownで後代がimputeされている家系の両親を補完する
# 後代は1サンプルのみ

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from ParentsImputerByProgeny import *


#################### VCFParentsKnownProgenyImputed ####################

class VCFParentsKnownProgenyImputed(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps_mat: list[list[int]],
					ref_haps_pat: list[list[int]],
					map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		self.imputer = ParentsImputerByProgeny(records, ref_haps_mat,
													ref_haps_pat, map_, 0.01)
		self.records: list[VCFFamilyRecord] = records
	
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
		self.imputer.impute()
