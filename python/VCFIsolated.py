from __future__ import annotations

# coding: utf-8
# VCFIsolated.py
# isolated samples

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFImputable import *
from Genotype import Genotype
from OptionSmall import *


#################### VCFIsolated ####################

# imputed samples at the beginning, followed by reference samples
class VCFIsolated(VCFImputable):
	def __init__(self, samples: list[str], num_imputed_samples: int,
						records: list[GenoRecord], map_: Map, vcf: VCFSmall):
		VCFImputable.__init__(self, samples, map_, vcf)
		self.records: list[GenoRecord] = records
		self.num_imputed_samples = num_imputed_samples
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	# divide VCF by 1 cM
	def divide_by_cM(self) -> Iterator[VCFIsolated]:
		# 1cMを超えても10個以内で10cM以内なら塊とみなす
		rs = self.records[:1]
		for i in range(1, len(self)):
			record = self.records[i]
			if self.is_block(record, rs):
				rs.append(record)
			else:
				yield VCFIsolated(self.samples, self.num_imputed_samples,
														rs, self.map, self.vcf)
				rs = [record]
		
		yield VCFIsolated(self.samples, self.num_imputed_samples,
														rs, self.map, self.vcf)
	
	def collect_haplotype_from_refs(self) -> list[Haplotype]:
		return [ self.clip_haplotype(i, j)
					for i in range(self.num_imputed_samples, self.num_samples())
					for j in (0, 1) ]
	
	def collect_haplotypes_mat(self, sample_index: int) -> list[Haplotype]:
		return self.collect_haplotype_from_refs()
	
	def collect_haplotypes_pat(self, sample_index: int) -> list[Haplotype]:
		return self.collect_haplotype_from_refs()
	
	def impute_cM(self, prev_haps: list[HaplotypePair]) -> list[HaplotypePair]:
		haps: list[HaplotypePair] = []
		for sample_index, prev_hap in enumerate(prev_haps):
			hap = self.impute_cM_each_sample(prev_hap, sample_index, True)
			haps.append(hap)
		return haps
	
	def impute(self) -> None:
		h = Haplotype.default_value()
		haps = [(h, h)] * self.num_imputed_samples
		for vcf_cM in self.divide_by_cM():
			haps = vcf_cM.impute_cM(haps)
		print("%d samples are imputed." % self.num_imputed_samples)
	
	@staticmethod
	def create(orig_vcf: VCFSmall, merged_vcf: VCFGeno,
						samples: list[str], references: list[str],
						op: OptionSmall) -> VCFIsolated:
		samples.sort()	# C++と合わせるため
		columns = orig_vcf.extract_columns(samples)
		ref_columns = merged_vcf.extract_columns(references)
		new_samples = [ orig_vcf.samples[c-9] for c in columns ] + references
		header = orig_vcf.trim_header(new_samples)
		new_records: list[GenoRecord] = []
		for i in range(len(orig_vcf)):
			record = orig_vcf.records[i]
			imputed_record = merged_vcf.records[i]
			pos = imputed_record.pos
			geno = ([ Genotype.all_gt_to_int(record.v[c]) for c in columns ] +
								[ imputed_record.geno[c] for c in ref_columns ])
			new_record = GenoRecord(pos, geno)
			new_records.append(new_record)
		vcf = VCFIsolated(new_samples, len(samples),
								new_records, op.map, orig_vcf)
		return vcf

__all__ = ['VCFIsolated']
