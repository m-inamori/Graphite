from __future__ import annotations

# coding: utf-8
# VCFImputable.py

from abc import ABC, abstractmethod
from itertools import *
from VCF import *
from Haplotype import *
from Map import *
from typing import Iterator, Tuple, TypeVar

R = TypeVar('R', bound='VCFRecord', covariant=True)


#################### VCFImputable ####################

class VCFImputable(VCFSmallBase, VCFMeasurable, ABC):
	def __init__(self, map_: Map):
		VCFSmallBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
	
	@abstractmethod
	def collect_haplotypes_mat(self, sample_index: int) -> list[Haplotype]:
		pass
	
	@abstractmethod
	def collect_haplotypes_pat(self, sample_index: int) -> list[Haplotype]:
		pass
	
	def get_int_gts(self, sample_index: int) -> list[int]:
		int_gts: list[int] = [ self.get_record(i).get_int_gt(sample_index)
												for i in range(len(self)) ]
		return int_gts
	
	def clip_haplotype(self, sample_index: int, side: int) -> Haplotype:
		hap = self.clip_raw_haplotype(sample_index, side)
		return Haplotype(hap, sample_index, side)
	
	def is_block(self, record: VCFRecord, rs: list[R]) -> bool:
		length = self.cM(record.pos()) - self.cM(rs[0].pos())
		if length < 1.0:
			return True
		elif len(rs) < 10 and length < 10.0:
			return True
		else:
			return False
	
	def impute_cM_each_sample(self, prev_hap: HaplotypePair,
								sample_index: int, exec: bool) -> HaplotypePair:
		int_gts: list[int] = self.get_int_gts(sample_index)
		haps_mat: list[Haplotype] = self.collect_haplotypes_mat(sample_index)
		haps_pat: list[Haplotype] = self.collect_haplotypes_pat(sample_index)
		seed: int = self.get_record(0).pos()
		hap = Haplotype.impute(int_gts, haps_mat, haps_pat, prev_hap, seed)
		if exec:	# actually impute?
			self.set_haplotype(hap, sample_index)
		return hap
	
	def set_haplotype(self, hap: HaplotypePair, sample_index: int):
		hap_mat, hap_pat = hap
		for i in range(len(self)):
			record = self.get_record(i)
			gt = "%d|%d" % (hap_mat.hap[i], hap_pat.hap[i])
			record.set_GT(sample_index, gt)
