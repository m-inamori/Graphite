
# coding: utf-8
# VCFSelfFillableRecord.py

from __future__ import annotations
from functools import reduce
from itertools import *
from collections import Counter

import VCF
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from TypeDeterminer import ParentComb
from common import *
from Genotype import Genotype


#################### VCFSelfFillableRecord ####################

class VCFSelfFillableRecord(GenoRecord):
	def __init__(self, pos: int, geno: list[int], index: int,
									type: SelfFillType, pair: ParentComb,
									probs: list[VCF.Probs]) -> None:
		super().__init__(pos, geno)
		self.index: int = index
		self.type: SelfFillType = type
		self.pair: ParentComb = pair
		self.probs: list[VCF.Probs] = probs
	
	def from_which_chrom(self, i: int, mat: bool) -> int:
		j = 0 if mat else 1
		parent_gt = self.geno[0]
		return 1 if (parent_gt & 1) == ((self.geno[i] >> j) & 1) else 2
	
	@staticmethod
	def convert(record: VCFImpSelfRecord,
						probs: list[VCF.Probs]) -> VCFSelfFillableRecord:
		type = record.get_fill_type()
		return VCFSelfFillableRecord(record.pos, record.geno, record.index,
													type, record.pair, probs)
