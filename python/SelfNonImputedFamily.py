from __future__ import annotations

# coding: utf-8
# SelfNonImputedFamily.py
# 自殖で一つもimputedが無い家系

from typing import Optional

from VCFFamily import *
from Map import *
from KnownFamily import *
from VCFSelfNoImputed import VCFSelfNoImputed
from VCFSelfNoImputedRough import VCFSelfNoImputedRough
from OptionSmall import OptionSmall


# HMMにrefを使っても計算量が十分に小さいか
def is_small(family: Family, ref_haps: list[list[int]],
									L: int, op: OptionSmall) -> bool:
	N = family.num_progenies()
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R: float = NH**2 * 4**N * (2*NH + 2*N - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def is_small_ref(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFSmall, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFSmallBase]:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		samples = [family.mat] + family.progenies
		vcf = VCFSmall.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		if is_small(family, ref_haps, len(families), op):
			vcf1 = VCFSelfNoImputed(vcf.header, vcf.records, ref_haps, op.map)
			vcf1.impute()
			vcfs.append(vcf1)
		elif is_small_ref(ref_haps, len(families), op):
			vcf2 = VCFSelfNoImputedRough(vcf.header, vcf.records,
														ref_haps, op.map)
			vcf2.impute()
			vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d self families have been imputed." % len(families))
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
