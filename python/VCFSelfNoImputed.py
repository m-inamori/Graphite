from __future__ import annotations

# coding: utf-8
# VCFSelfNoImputed.py
# 自殖でimputedなサンプルが一つもない

from VCF import VCFSmall
from VCFSelfImputable import VCFSelfImputable
from GenoRecord import GenoRecord
from SelfImputer import SelfImputer
from Map import *


#################### VCFSelfNoImputed ####################

class VCFSelfNoImputed(VCFSelfImputable):
	def __init__(self, samples: list[str],
						records: list[GenoRecord],
						ref_haps: list[list[int]],
						map_: Map, vcf: VCFSmall) -> None:
		VCFSelfImputable.__init__(self, samples, vcf)
		self.records: list[GenoRecord] = records
		self.imputer = SelfImputer(records, ref_haps, map_, 0.01)
	
	##### virtual methods for VCFGenoBase #####
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_records(self) -> list[GenoRecord]:
		return [ r for r in self.records ]
	
	##### virtual methods for VCFGenoBase #####
	def impute(self) -> None:
		self.imputer.impute()
