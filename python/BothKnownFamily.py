from __future__ import annotations

# coding: utf-8
# BothKnownFamily.py
# 両親が分っているがphasingされていない家系を補完する

from functools import reduce
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from Map import *
from VCFBothKnown import VCFBothKnown
from OptionSmall import OptionSmall


def is_small(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
			families: list[Family], op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		if is_small(ref_haps, len(families), op):
			vcf1 = VCFBothKnown(vcf.samples, vcf.records,
										ref_haps, op.map, orig_vcf)
			vcf1.impute()
			vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose parents are known have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
