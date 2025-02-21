from __future__ import annotations

# coding: utf-8
# VCFHeteroHomoOnePhased.py
# 片親がphasingされて片親は分っているがphasingされていないヘテロ×ホモファミリー

from abc import abstractmethod
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from VCFFillableRecord import *
from Map import *


#################### VCFHeteroHomoOnePhased ####################

class VCFHeteroHomoOnePhased(VCFBase, VCFSmallBase,
								VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]],
						records: list[VCFFillableRecord],
						is_mat_hetero: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFillableRecord] = records
		self.is_mat_hetero = is_mat_hetero
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	@abstractmethod
	def impute(self) -> None:
		pass
