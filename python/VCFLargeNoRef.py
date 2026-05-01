from __future__ import annotations

# coding: utf-8
# VCFLargeNoRef.py
# 両親ともにRefになく後代ともphasedとN/A

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFImputable import VCFImputable
from ParentNoRefImputer import *
from ParentRefImputer import *
from ProgenyRefImputer import *
from Map import *
from Genotype import Genotype


#################### VCFLargeNoRef ####################

class VCFLargeNoRef(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
													map_: Map, ref_vcf: VCFGeno):
		VCFImputable.__init__(self, samples, ref_vcf.vcf)
		self.records: list[VCFFamilyRecord] = records
		self.mat_imputer = ParentNoRefImputer(records, map_,
												True, ref_vcf, 0.01)
		self.pat_imputer = ParentRefImputer(records, map_, True, ref_vcf, 0.01)
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
		self.mat_imputer.impute()
		self.pat_imputer.impute()
		for i in range(self.num_progenies()):
			self.prog_imputer.impute(i)
