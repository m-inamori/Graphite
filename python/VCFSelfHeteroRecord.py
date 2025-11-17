from __future__ import annotations

# coding: utf-8
# VCFSelfHeteroRecord.py
# 自殖で親がヘテロ

from collections import defaultdict, Counter
from math import log10
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from VCF import VCFRecord, VCFBase, VCFSmallBase
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from TypeDeterminer import ParentComb
from Map import *
import Imputer
from Genotype import Genotype
from option import *
from inverse_graph import *
from graph import Node
from invgraph import InvGraph


#################### VCFSelfHeteroRecord ####################

class VCFSelfHeteroRecord(VCFImpSelfRecord):
	def __init__(self, pos: int, geno: list[int],
							i: int, parents_wrong_type: str,
							pair: ParentComb) -> None:
		super().__init__(pos, geno, i, parents_wrong_type, pair)
	
	def num_progenies(self) -> int:
		return len(self.geno) - 1
	
	def is_imputable(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def get_fill_type(self) -> SelfFillType:
		if self.is_imputable():
			return SelfFillType.P01
		else:
			return SelfFillType.IMPUTABLE
	
	def set_haplo(self, h: int) -> None:
		self.geno[0] = Genotype.PH_01 if h == 0 else Genotype.PH_10
