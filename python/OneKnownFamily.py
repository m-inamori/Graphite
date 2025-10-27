from __future__ import annotations

# coding: utf-8
# OneKnownFamily.py
# 片親がphasingされていないがknownでもう片親はunknownな家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from Map import *
from KnownFamily import *
from VCFOneParentKnown import VCFOneParentKnown
from OptionSmall import OptionSmall


# HMMにrefを使っても計算量が十分に小さいか
def is_small(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		if is_small(ref_haps, len(families), op):
			vcf1 = VCFOneParentKnown(vcf.samples, vcf.records, ref_haps,
											family.mat_known, op.map, orig_vcf)
			vcf1.impute_known_parent()
			parent = vcf1.samples[0] if family.mat_known else vcf1.samples[1]
			vcf2 = vcf1.extract_by_samples([parent])
			vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is known and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
