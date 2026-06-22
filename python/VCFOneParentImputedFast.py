from __future__ import annotations

# coding: utf-8
# VCFOneParentImputedFast.py
# 片側だけimputeされている
# ヘテロ×ホモを手掛かりにimputeする

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFFillable import VCFFillable
from VCFSmallFillable import VCFSmallFillable
from VCFImputable import VCFImputable
from VCFImpHeteroHomo import VCFImpHeteroHomo
from VCFHeteroImpHomo import VCFHeteroImpHomo
from VCFFillableRecord import VCFFillableRecord
from ParentImputer import *
from ProgenyImputer import *
from ClassifyRecord import FillType, classify_family_record
from Map import *


#################### VCFOneParentImputedFast ####################

class VCFOneParentImputedFast(VCFImputable):
	def __init__(self, samples: list[str], records: list[VCFFamilyRecord],
					should_impute_mat: bool, map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, vcf)
		self.records: list[VCFFamilyRecord] = records
		self.should_impute_mat: bool = should_impute_mat
		self.gmap = map_
	
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
	def classify_records(self) -> list[list[VCFFillableRecord]]:
		cols = self.vcf.extract_columns(self.samples)
		# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
		rss: list[list[VCFFillableRecord]] = [ [] for _ in range(4) ]
		for index, record in enumerate(self.records):
			pair, type = classify_family_record(record)
			ref_record = self.vcf.records[index]
			probs = ref_record.parse_PL(record.geno, cols)
			rss[type.value].append(VCFFillableRecord(record.pos, record.geno,
														index, type, pair, probs))
		return rss
	
	def create(self, rs: list[VCFFillableRecord],
								is_mat_hetero: bool) -> VCFImputable:
		if is_mat_hetero != self.should_impute_mat:
			return VCFImpHeteroHomo(self.samples, rs, is_mat_hetero,
														self.gmap, self.vcf)
		else:
			return VCFHeteroImpHomo(self.samples, rs, is_mat_hetero,
														self.gmap, self.vcf)
	
	def merge_vcf(self, rss: list[list[VCFFillableRecord]]) -> VCFFillable:
		rs = rss[0] + rss[1] + rss[2] + rss[3]
		rs.sort(key=lambda r: r.pos)
		return VCFSmallFillable(self.samples, rs, self.vcf)
	
	##### virtual methods for VCFImputable #####
	def impute(self) -> None:
		rss = self.classify_records()
		mat_vcf = self.create(rss[FillType.MAT.value], True)
		mat_vcf.impute()
		pat_vcf = self.create(rss[FillType.PAT.value], False)
		pat_vcf.impute()
		merged_vcf = self.merge_vcf(rss)
		merged_vcf.modify(False)
