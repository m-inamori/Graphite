from __future__ import annotations

# coding: utf-8
# VCFHeteroHomoOnePhased.py
# 片親がphasingされて片親は分っているがphasingされていないヘテロ×ホモファミリー

from abc import abstractmethod
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFFillableRecord import *
from Map import *


#################### VCFHeteroHomoOnePhased ####################

class VCFHeteroHomoOnePhased(VCFFamilyBase, VCFMeasurable):
	def __init__(self, samples: list[str],
						records: list[VCFFillableRecord],
						is_mat_hetero: bool, map_: Map, vcf: VCFSmall):
		VCFFamilyBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.samples: list[str] = samples
		self.records: list[VCFFillableRecord] = records
		self.is_mat_hetero = is_mat_hetero
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	@abstractmethod
	def impute(self) -> None:
		pass
