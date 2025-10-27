from __future__ import annotations

# coding: utf-8
# VCFFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCF import *
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFSelfHetero import VCFSelfHetero
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from VCFSelfFillableRecord import VCFSelfFillableRecord
from VCFHeteroHomo import *
from SelfGroups import SelfGroups
from SelfRecordSet import SelfRecordSet
from Genotype import Genotype
from option import *
from common import *


#################### VCFSelfFillable ####################

class VCFSelfFillable(VCFGenoBase):
	def __init__(self, samples: list[str],
						records: list[VCFSelfFillableRecord],
						vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
		self.records: list[VCFSelfFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def modify(self) -> None:
		# FillTypeでrecordを分ける
		groups = SelfGroups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.determine_parents_phasing()
	
	@staticmethod
	def fill(vcfs: list[VCFSelfHetero],
					records: list[VCFImpSelfRecord]) -> VCFSelfFillable:
		merged_records = VCFSelfFillable.merge_records(vcfs, records)
		vcf = VCFSelfFillable(vcfs[0].samples, merged_records, vcfs[0].vcf)
		vcf.modify()
		return vcf
	
	@staticmethod
	def merge_records(vcfs: list[VCFSelfHetero], records: list[VCFImpSelfRecord]
											) -> list[VCFSelfFillableRecord]:
		all_records = [ record for record in chain(
							(r for vcf in vcfs for r in vcf.records), records) ]
		all_records.sort(key=lambda record: record.index)
		probs = VCFSelfFillable.calc_probs(all_records, vcfs[0])
		new_records = [ VCFSelfFillableRecord.convert(r, p)
									for r, p in zip(all_records, probs) ]
		return new_records
	
	@staticmethod
	def calc_probs(records: list[VCFImpSelfRecord],
					vcf: VCFGenoBase) -> list[list[tuple[float, float, float]]]:
		orig_vcf = vcf.vcf
		cols = orig_vcf.extract_columns(vcf.samples)
		prob_table = []
		for record, orig_record in zip(records, orig_vcf.records):
			probs = orig_record.parse_PL(record.geno, cols)
			prob_table.append(probs)
		return prob_table
