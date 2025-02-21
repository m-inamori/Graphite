from __future__ import annotations

# coding: utf-8
# OnePhasedFamily.py
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

from functools import reduce
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from Map import *
from VCFNoParentImputed import VCFNoParentImputed
from VCFOneParentImputedRough import VCFOneParentImputedRough

def is_small(ref_haps: list[list[int]], L: int) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1)
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFSmall, ref_haps: list[list[int]],
							families: list[Family], gmap: Map) -> VCFSmallBase:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		if is_small(ref_haps, len(families)):
			vcf1 = VCFNoParentImputed(vcf.header, vcf.records, ref_haps, gmap)
			vcf1.impute()
			vcfs.append(vcf1)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
