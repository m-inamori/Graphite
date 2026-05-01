from __future__ import annotations

# coding: utf-8
# VCFImputable.py

from abc import abstractmethod

from VCF import VCFSmall
from VCFGeno import VCFGenoBase


#################### VCFSelfImputable ####################

class VCFSelfImputable(VCFGenoBase):
	def __init__(self, samples: list[str], vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
	
	def num_progenies(self) -> int:
		return self.num_samples() - 1
	
	@abstractmethod
	def impute(self) -> None:
		pass
