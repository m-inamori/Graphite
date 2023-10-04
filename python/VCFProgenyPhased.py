from __future__ import annotations

# coding: utf-8
# VCFProgenyPhased.py
# 子どもが一つでもphasingされた家系
# 親だけphasingする

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from VCFImputable import *
from Map import *
import Imputer


#################### VCFProgenyPhased ####################

class VCFProgenyPhased(VCFBase, VCFFamilyBase, VCFImputable):
	def __init__(self, header: list[list[str]],
							records: list[VCFFamilyRecord],
							map_: Map, ref_vcf: VCFSmall):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		VCFImputable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_vcf: VCFSmall = ref_vcf
		self.selection = 0
	
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
	def divide_by_cM(self) -> Iterator[VCFProgenyPhased]:
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
				yield VCFProgenyPhased(self.header, rs, self.map, ref_vcf_cM)
				rs = [record]
				ref_rs = [ref_record]
		
		ref_vcf_cM = VCFSmall(self.ref_vcf.header, ref_rs)
		yield VCFProgenyPhased(self.header, rs, self.map, ref_vcf_cM)
	
	def collect_haplotypes_from_phased_progeny(self, j: int) -> list[Haplotype]:
		return [ self.clip_haplotype(2, j) ]
	
	def collect_haplotype_from_refs(self) -> list[Haplotype]:
		return [ self.clip_ref_haplotype(i, j)
					for i in range(self.ref_vcf.num_samples())
					for j in (0, 1) ]
	
	def collect_haplotypes_mat(self, sample_index: int) -> list[Haplotype]:
		i = (self.selection + sample_index) % 2
		return self.collect_haplotypes_from_phased_progeny(i)
	
	def collect_haplotypes_pat(self, sample_index: int) -> list[Haplotype]:
		return self.collect_haplotype_from_refs()
	
	def impute_cM(self, prev_haps: list[HaplotypePair],
							sample_indices: list[int]) -> list[HaplotypePair]:
		haps: list[HaplotypePair] = []
		for sample_index, prev_hap in zip(sample_indices, prev_haps):
			hap = self.impute_cM_each_sample(prev_hap, sample_index, False)
			haps.append(hap)
		return haps
	
	def impute_core(self, vcf_cMs: list[VCFProgenyPhased]
											) -> list[list[HaplotypePair]]:
		sample_indices = [ i for i, s in enumerate(vcf_cMs[0].samples[:2])
																if s != '0' ]
		N = len(sample_indices)
		whole_haplotype_pair: list[list[HaplotypePair]] = [[] for _ in range(N)]
		h = Haplotype.default_value()
		haps = [(h, h)] * len(sample_indices)
		for vcf_cM in vcf_cMs:
			haps = vcf_cM.impute_cM(haps, sample_indices)
			for i in range(N):
				whole_haplotype_pair[i].append(haps[i])
		return whole_haplotype_pair
	
	def score_whole(self, haps: list[HaplotypePair],
					sample_index: int, vcf_cMs: list[VCFProgenyPhased]) -> int:
		score = 0
		int_gts = self.get_int_gts(sample_index)
		first = 0
		for i in range(len(vcf_cMs)):
			vcf_cM = vcf_cMs[i]
			hap = haps[i]
			last = first + len(vcf_cM)
			score += Haplotype.score(int_gts[first:last], hap[0], hap[1])
			first = last
		return score
	
	def inverse_selection(self):
		self.selection = 1 if self.selection == 0 else 0
	
	def impute(self):
		vcf_cMs = list(self.divide_by_cM())
		hap_pairs1 = self.impute_core(vcf_cMs)
		for vcf_cM in vcf_cMs:
			vcf_cM.inverse_selection()
		hap_pairs2 = self.impute_core(vcf_cMs)
		sample_indices = [ i for i, s in enumerate(vcf_cMs[0].samples[:2])
																if s != '0' ]
		score1 = sum(self.score_whole(hap_pairs1[i], i, vcf_cMs)
												for i in sample_indices)
		score2 = sum(self.score_whole(hap_pairs2[i], i, vcf_cMs)
												for i in sample_indices)
		_, hap_pairs = max([(score1, hap_pairs1), (score2, hap_pairs2)],
														key=lambda v: v[0])
		for j, vcf_cM in enumerate(vcf_cMs):
			for i in sample_indices:
				vcf_cM.set_haplotype(hap_pairs[i][j], i)
	
	@staticmethod
	def create(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
							family_samples: list[str], phased_index: int,
							map_: Map, ref_vcf: VCFSmall) -> VCFProgenyPhased:
		# extract parents and the phased progeny
		samples = family_samples[:2] + [family_samples[phased_index]]
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf, orig_vcf, samples)
		return VCFProgenyPhased(vcf.header, vcf.records, map_, ref_vcf)

__all__ = ['VCFProgenyPhased']
