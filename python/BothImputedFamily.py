from __future__ import annotations

# coding: utf-8
# BothImputedFamily.py
# 両親がimputeされている家系を補完する

from VCFFamily import *
from VCFBothParentImputed import *
from Map import *
from KnownFamily import KnownFamily
from typing import Optional

def impute(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
				families: list[KnownFamily], gmap: Map) -> Optional[VCFSmall]:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf, orig_vcf,
														family.samples())
		family_vcf = VCFBothParentImputed(vcf.header, vcf.records, gmap)
		family_vcf.impute()
		vcfs.append(family_vcf)
	
	if not vcfs:
		return None
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	print("%d both parent imputed families have been imputed." % len(vcfs))
	return new_vcf

__all__ = ['impute']
