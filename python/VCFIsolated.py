from __future__ import annotations

# coding: utf-8
# VCFIsolated.py
# isolated samples

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from Haplotype import *
import RefCommon
from Genotype import Genotype
from OptionSmall import *


#################### VCFIsolated ####################

# imputed samples at the beginning, followed by reference samples
class VCFIsolated(VCFGenoBase, VCFMeasurable):
	def __init__(self, samples: list[str], num_imputed_samples: int,
						records: list[GenoRecord], map_: Map, vcf: VCFSmall):
		VCFGenoBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records: list[GenoRecord] = records
		self.num_imputed_samples = num_imputed_samples
	
	##### virtual methods for VCFGenoBase #####
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_records(self) -> list[GenoRecord]:
		return [ r for r in self.records ]
	
	##### non-virtual methods #####
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
	
	def clip_haplotype(self, sample_index: int, side: int) -> Haplotype:
		hap = self.clip_raw_haplotype(sample_index, side)
		return Haplotype(hap, sample_index, side)
	
	def is_block(self, record: GenoRecord, rs: list[GenoRecord]) -> bool:
		length = self.cM(record.pos) - self.cM(rs[0].pos)
		if length < 1.0:
			return True
		elif len(rs) < 10 and length < 10.0:
			return True
		else:
			return False
	
	def get_int_gts(self, sample_index: int) -> list[int]:
		int_gts: list[int] = [ self.get_record(i).unphased(sample_index)
												for i in range(len(self)) ]
		return int_gts
	
	def impute_cM_each_sample(self, prev_hap: HaplotypePair,
								sample_index: int, exec: bool) -> HaplotypePair:
		int_gts: list[int] = self.get_int_gts(sample_index)
		haps_mat: list[Haplotype] = self.collect_haplotypes_mat(sample_index)
		haps_pat: list[Haplotype] = self.collect_haplotypes_pat(sample_index)
		seed: int = self.get_record(0).pos
		hap = Haplotype.impute(int_gts, haps_mat, haps_pat, prev_hap, seed)
		if exec:	# actually impute?
			self.set_haplotype(hap, sample_index)
		return hap
	
	def set_haplotype(self, hap: HaplotypePair, sample_index: int) -> None:
		hap_mat, hap_pat = hap
		for i in range(len(self)):
			record = self.get_record(i)
			gt = hap_mat.hap[i] | (hap_pat.hap[i] << 1) | 4
			record.geno[sample_index] = gt
	
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
	
	@staticmethod
	def create_ref(orig_vcf: VCFSmall, phased_vcf: VCFGeno,
						samples: list[str], references: list[str],
						op: OptionSmall) -> VCFIsolated:
		vcf = VCFGeno.extract_samples(samples, orig_vcf);
		new_samples = samples + references
		records = RefCommon.merge_records(phased_vcf, vcf, new_samples)
		return VCFIsolated(new_samples, len(samples),
									records, op.map, orig_vcf)

__all__ = ['VCFIsolated']
