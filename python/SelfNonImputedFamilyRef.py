from __future__ import annotations

# coding: utf-8
# SelfFamily.py
# 自殖の家系

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
import SelfNonImputedFamily
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


#################### SelfNonImputedFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno, ref_haps: list[list[int]],
										families: list[KnownFamily],
										op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_records(phased_vcf, vcf, family.samples())
		vcf1 = SelfNonImputedFamily.create_family_vcf(family, records,
														len(families), ref_haps,
														orig_vcf, op)
		if vcf1 is not None:
			vcf1.impute()
			vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	
	print("%d self families have been imputed." % len(families))
	return new_vcf

__all__ = ['impute']
