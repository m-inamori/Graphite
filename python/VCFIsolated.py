from __future__ import annotations

# coding: utf-8
# VCFIsolated.py
# isolated samples

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCF import *
from VCFImputable import *
from OptionSmall import *


#################### VCFIsolated ####################

# imputed samples at the beginning, followed by reference samples
class VCFIsolated(VCFBase, VCFImputable):
	def __init__(self, header: list[list[str]], num_imputed_samples: int,
									records: list[VCFRecord], map_: Map):
		VCFBase.__init__(self, header)
		VCFImputable.__init__(self, map_)
		self.records: list[VCFRecord] = records
		self.num_imputed_samples = num_imputed_samples
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
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
				yield VCFIsolated(self.header, self.num_imputed_samples,
															rs, self.map)
				rs = [record]
		
		yield VCFIsolated(self.header, self.num_imputed_samples,
															rs, self.map)
	
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
	
	def extract_isolated_samples(self) -> VCFSmall:
		isolated_samples = self.samples[:self.num_imputed_samples]
		return self.extract_samples(isolated_samples)
	
	@staticmethod
	def create(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
						samples: list[str], references: list[str],
						op: OptionSmall) -> VCFIsolated:
		sample_columns = orig_vcf.extract_columns(samples)
		ref_columns = merged_vcf.extract_columns(references)
		new_samples = [ orig_vcf.samples[c-9]
								for c in sample_columns ] + references
		header = orig_vcf.trim_header(new_samples)
		new_records = []
		for i in range(len(orig_vcf)):
			record = orig_vcf.records[i]
			imputed_record = merged_vcf.records[i]
			v = (record.v[:9] + [ record.v[c] for c in sample_columns ] +
								[ imputed_record.v[c] for c in ref_columns ])
			new_record = VCFRecord(v, new_samples)
			new_records.append(new_record)
		vcf = VCFIsolated(header, len(samples), new_records, op.map)
		return vcf

__all__ = ['VCFIsolated']
