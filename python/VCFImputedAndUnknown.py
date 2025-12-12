from __future__ import annotations

# coding: utf-8
# VCFImputedAndUnknown.py
# 片親がimputeされていて、もう片親はUnknown
# まず片親

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from ProgenyImputerByOneParent import *
from ParentImputerByProgeny import *
from ProgenyImputer import *
from Map import *


#################### VCFImputedAndUnknown ####################

class VCFImputedAndUnknown(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool,
					map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.prog_imputer = ProgenyImputerByOneParent(records, ref_haps,
														is_mat_imputed,
														map_, 0.01)
		self.parent_imputer = ParentImputerByProgeny(records, ref_haps,
														is_mat_imputed,
														map_, 0.01)
		self.progs_imputer = ProgenyImputer(records, map_, 0.01)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def impute(self) -> None:
		# 最初の後代をimpute
		self.prog_imputer.impute(0)
		# 親をimpute
		self.parent_imputer.impute()
		# 残りの後代をimpute
		for i in range(1, self.num_progenies()):
			self.progs_imputer.impute(i)
