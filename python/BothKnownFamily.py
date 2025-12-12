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
from VCFBothKnown import VCFBothKnown
from Map import *
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall


def is_small(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

# is_smallが通るギリギリのNHを求める
def compute_upper_NH(family: Family, M: int, L: int, op: OptionSmall) -> int:
	for NH in count(2):
		R = NH**2 * (2*NH - 1) / op.precision_ratio
		if not (R * M < 10**8 and R < 10**5 and L * R * M < 10**9):
			break
	return NH - 1

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
			families: list[Family], op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		NH = compute_upper_NH(family, len(vcf), len(families), op)
		if is_small(ref_haps, len(families), op):
			vcf1 = VCFBothKnown(vcf.samples, vcf.records,
								ref_haps, ref_haps, op.map, orig_vcf)
			vcf1.impute()
			vcfs.append(vcf1)
		elif NH >= 10:
			NH2 = min(20, NH)
			gts_mat = [ r.geno[0] for r in vcf.records ]
			gts_pat = [ r.geno[1] for r in vcf.records ]
			ref_haps_mat = filter_haplotypes(ref_haps, gts_mat, NH2)
			ref_haps_pat = filter_haplotypes(ref_haps, gts_pat, NH2)
			vcf2 = VCFBothKnown(vcf.samples, vcf.records,
								ref_haps_mat, ref_haps_pat, op.map, orig_vcf)
			vcf2.impute()
			vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d families whose parents are known have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
