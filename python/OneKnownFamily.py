from __future__ import annotations

# coding: utf-8
# OneKnownFamily.py
# 片親がphasingされていないがknownでもう片親はunknownな家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFImputable import VCFImputable
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

def create_family_vcf(family: KnownFamily, records: list[VCFFamilyRecord],
					num_families: int, ref_haps: list[list[int]],
					vcf: VCFSmall, op: OptionSmall) -> Optional[VCFImputable]:
	NH = compute_upper_NH(family, len(vcf), num_families, op)
	if is_small(ref_haps, num_families, op):
		return VCFOneParentKnown(family.samples(), records, ref_haps,
											family.mat_known, op.map, vcf)
	elif NH >= 10:
		NH2 = min(20, NH)
		if family.mat_known:
			gts = [ r.geno[0] for r in records ]
			ref_haps = filter_haplotypes(ref_haps, gts, NH2)
		else:
			gts = [ r.geno[1] for r in records ]
			ref_haps = filter_haplotypes(ref_haps, gts, NH2)
		return VCFOneParentKnown(family.samples(), records, ref_haps,
											family.mat_known, op.map, vcf)
	else:
		return None

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		parent = vcf.mat() if family.mat_known else vcf.pat()
		vcf1 = create_family_vcf(family, vcf.records, len(families),
													ref_haps, orig_vcf, op)
		if vcf1 is None:
			continue
		vcf1.impute()
		vcf1_parent = vcf1.extract_by_samples([parent])
		vcfs.append(vcf1_parent)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is known and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
