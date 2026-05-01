from __future__ import annotations

# coding: utf-8
# VCFImputable.py

from abc import abstractmethod

from VCF import VCFSmall
from VCFFamily import VCFFamilyBase


#################### VCFImputable ####################

class VCFImputable(VCFFamilyBase):
	def __init__(self, samples: list[str], vcf: VCFSmall) -> None:
		VCFFamilyBase.__init__(self, samples, vcf)
	
	@abstractmethod
	def impute(self) -> None:
		pass
