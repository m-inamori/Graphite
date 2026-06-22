from __future__ import annotations

# coding: utf-8
# VCFOneParentProgenyImputed.py
# 片親だけimputeされていてもう片親はknown

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFFamily import VCFFamilyRecord
from VCFImputable import VCFImputable
from ImputerByParentProgeny import ImputerByParentProgeny
from Map import *


#################### VCFOneParentProgenyImputed ####################

class VCFOneParentProgenyImputed(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], should_impute_mat: bool,
					map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.imputer = ImputerByParentProgeny(records, ref_haps,
												should_impute_mat, map_, 0.01)
		self.should_impute_mat = should_impute_mat
	
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
		self.imputer.impute()
