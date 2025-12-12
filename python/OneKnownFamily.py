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
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall


# HMMにrefを使っても計算量が十分に小さいか
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
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		NH = compute_upper_NH(family, len(vcf), len(families), op)
		parent = vcf.samples[0] if family.mat_known else vcf.samples[1]
		if is_small(ref_haps, len(families), op):
			vcf1 = VCFOneParentKnown(vcf.samples, vcf.records,
										ref_haps, ref_haps,
										family.mat_known, op.map, orig_vcf)
			vcf1.impute_known_parent()
			vcf1_parent = vcf1.extract_by_samples([parent])
			vcfs.append(vcf1_parent)
		elif NH >= 10:
			NH2 = min(20, NH)
			gts_mat = [ r.geno[0] for r in vcf.records ]
			gts_pat = [ r.geno[1] for r in vcf.records ]
			ref_haps_mat = filter_haplotypes(ref_haps, gts_mat, NH2)
			ref_haps_pat = filter_haplotypes(ref_haps, gts_pat, NH2)
			ref_haps1 = ref_haps_mat if family.mat_known else ref_haps_pat
			ref_haps2 = ref_haps_pat if family.mat_known else ref_haps_mat
			vcf2 = VCFOneParentKnown(vcf.samples, vcf.records,
										ref_haps1, ref_haps2,
										family.mat_known, op.map, orig_vcf)
			vcf2.impute_known_parent()
			vcf2_parent = vcf2.extract_by_samples([parent])
			vcfs.append(vcf2_parent)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is known and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
