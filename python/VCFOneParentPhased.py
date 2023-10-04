from __future__ import annotations

# coding: utf-8
# VCFOneParentPhased.py
# 片方の親だけがphasingされた家系

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from VCFImputable import *
from Map import *
import Imputer
from option import *


#################### VCFHOneParentPhased ####################

class VCFOneParentPhased(VCFBase, VCFFamilyBase, VCFImputable):
	def __init__(self, header: list[list[str]],
						records: list[VCFFamilyRecord],
						mat_p: bool, map_: Map, ref_vcf: VCFSmall):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		VCFImputable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.is_mat_phased: bool = mat_p
		self.ref_vcf: VCFSmall = ref_vcf
	
	def get_header(self) -> list[list[str]]:
		return self.header
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def clip_ref_haplotype(self, sample_id: int, side: int) -> Haplotype:
		hap = self.ref_vcf.clip_raw_haplotype(sample_id, side)
		return Haplotype(hap, sample_id, side)
	
	# divide VCF by 1 cM
	def divide_by_cM(self) -> Iterator[VCFOneParentPhased]:
		# 1cMを超えても10個以内で10cM以内なら塊とみなす
		rs = self.records[:1]
		ref_rs = self.ref_vcf.records[:1]
		for i in range(1, len(self)):
			record = self.records[i]
			ref_record = self.ref_vcf.records[i]
			if self.is_block(record, rs):
				rs.append(record)
				ref_rs.append(ref_record)
			else:
				ref_vcf_cM = VCFSmall(self.ref_vcf.header, ref_rs)
				yield VCFOneParentPhased(self.header, rs, self.is_mat_phased,
														self.map, ref_vcf_cM)
				rs = [record]
				ref_rs = [ref_record]
		
		ref_vcf_cM = VCFSmall(self.ref_vcf.header, ref_rs)
		yield VCFOneParentPhased(self.header, rs, self.is_mat_phased,
														self.map, ref_vcf_cM)
	
	def collect_haplotypes_from_parents(self) -> list[Haplotype]:
		i = 0 if self.is_mat_phased else 1
		return [ self.clip_haplotype(i, j) for j in (0, 1) ]
	
	def collect_haplotype_from_refs(self) -> list[Haplotype]:
		return [ self.clip_ref_haplotype(i, j)
					for i in range(self.ref_vcf.num_samples())
					for j in (0, 1) ]
	
	def collect_haplotypes_mat(self, sample_index: int) -> list[Haplotype]:
		if self.is_mat_phased:
			return self.collect_haplotypes_from_parents()
		else:
			return self.collect_haplotype_from_refs()
	
	def collect_haplotypes_pat(self, sample_index: int) -> list[Haplotype]:
		if self.is_mat_phased:
			return self.collect_haplotype_from_refs()
		else:
			return self.collect_haplotypes_from_parents()
	
	def impute_cM(self, prev_haps: list[HaplotypePair]) -> list[HaplotypePair]:
		haps: list[HaplotypePair] = []
		for sample_index, prev_hap in enumerate(prev_haps, 2):
			hap = self.impute_cM_each_sample(prev_hap, sample_index, True)
			haps.append(hap)
		return haps
	
	def impute(self):
		if not self.records:
			return
		
		h = Haplotype.default_value()
		haps = [(h, h)] * (self.num_progenies())
		for vcf_cM in self.divide_by_cM():
			haps = vcf_cM.impute_cM(haps)
	
	@staticmethod
	def create(samples: list[str], is_mat_phased: bool,
				imputed_vcf: VCFSmallBase, orig_vcf: VCFSmallBase,
				gmap: Map, ref_vcf: VCFSmall) -> VCFOneParentPhased:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		new_vcf = VCFOneParentPhased(vcf.header, vcf.records,
											is_mat_phased, gmap, ref_vcf)
		return new_vcf
	
	# samplesに対応するRecordを作る
	@staticmethod
	def merge_records(vcfs: list[VCFFamilyBase], i: int, samples: list[str],
					dic_pos: Dict[str, tuple[VCFFamilyBase, int]]) -> VCFRecord:
		v = vcfs[0].get_family_record(i).v[:9]
		for sample in samples:
			vcf, c = dic_pos[sample]
			v.append(vcf.get_family_record(i).v[c])
		return VCFRecord(v, samples)

__all__ = ['VCFOneParentPhased']
