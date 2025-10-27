from __future__ import annotations

# coding: utf-8
# VCFImpSelfRecord.py
# VCFSelfFillableでfillされるRecord

from abc import ABCMeta, abstractmethod
from collections import defaultdict
from operator import xor
from enum import Enum

from GenoRecord import GenoRecord
from TypeDeterminer import ParentComb
from common import is_all_same


#################### SelfFillType ####################

class SelfFillType(Enum):
	P01 = 0		# 親のGenotypeは0/1
	FILLED = 1
	IMPUTABLE = 2
	UNABLE = 3


#################### VCFImpSelfRecord ####################

class VCFImpSelfRecord(GenoRecord, metaclass=ABCMeta):
	def __init__(self, pos: int, geno: list[int],
						index: int, parents_wrong_type: str, pair: ParentComb):
		super().__init__(pos, geno)
		self.index: int = index
		self.parents_wrong_type: str = parents_wrong_type
		self.pair: ParentComb = pair
	
	def is_fixed(self) -> bool:
		return (self.pair in (ParentComb.P00x00, ParentComb.P11x11) or
				(self.pair.is_heterohomo() and self.is_right()))
	
	def is_right(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def enable_modification(self) -> None:
		self.parents_wrong_type = 'Modifiable'
	
	@abstractmethod
	def is_imputable(self) -> bool:
		pass
	
	@abstractmethod
	def get_fill_type(self) -> SelfFillType:
		pass
