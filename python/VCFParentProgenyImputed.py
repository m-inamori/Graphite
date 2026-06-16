from __future__ import annotations

# coding: utf-8
# VCFParentProgenyImputed.py
# 片親がimputeされていて、もう片親はUnknown
# まず片親

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImputable import VCFImputable
from ProgenyImputerByOneParent import *
from ParentImputerByProgeny import *
from ProgenyImputer import *
from Map import *


#################### VCFParentProgenyImputed ####################

class VCFParentProgenyImputed(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], should_impute_mat: bool,
					map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.imputer = ParentImputerByProgeny(records, ref_haps,
														should_impute_mat,
														map_, 0.01)
	
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
