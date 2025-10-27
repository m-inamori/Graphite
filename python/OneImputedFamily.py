from __future__ import annotations

# coding: utf-8
# OneImputedFamily.py
# 片親がimputeされていてもう片親はunknownな家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from Map import *
from VCFOneParentImputed import VCFOneParentImputed
from VCFOneParentImputedRough import VCFOneParentImputedRough
from KnownFamily import *
from OptionSmall import OptionSmall


# HMMにrefを使っても計算量が十分に小さいか
def is_small(family: Family, ref_haps: list[list[int]], L: int,
												op: OptionSmall) -> bool:
	N = family.num_progenies()
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R: int = NH**2 * 4**N * (2*NH + 2*N - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

# 後代はHMM
def is_small_ref(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		for prog in family.progenies:
			samples = [family.mat, family.pat, prog]
			vcf = VCFFamily.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
			if is_small(family, ref_haps, len(families), op):
				vcf1 = VCFOneParentImputed(samples, vcf.records, ref_haps,
											family.mat_known, op.map, orig_vcf)
				vcf1.impute()
				vcfs.append(vcf1)
			elif is_small_ref(ref_haps, len(families), op):
				vcf2 = VCFOneParentImputedRough(samples, vcf.records, ref_haps,
											family.mat_known, op.map, orig_vcf)
				vcf2.impute()
				vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is imputed and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
