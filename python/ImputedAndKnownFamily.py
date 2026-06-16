from __future__ import annotations

# coding: utf-8
# ImputedAndKnownFamily.py
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

from typing import Optional, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFImputable import VCFImputable
from VCFOneParentImputed import VCFOneParentImputed
from VCFOneParentImputedRough import VCFOneParentImputedRough
from VCFOneParentImputedFast import VCFOneParentImputedFast
from Map import *
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall


#################### ImputedAndKnownFamily ####################

# HMMにrefを使っても計算量が十分に小さいか
def is_small(family: Family, ref_haps: list[list[int]],
									L: int, op: OptionSmall) -> bool:
	N = family.num_progenies()
	if N > 2:
		return False
	
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R: float = NH**2 * 4**N * (2*NH + 2*N - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def is_small_ref(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

# is_small_refが通るギリギリのNHを求める
def compute_upper_NH(family: Family, M: int, L: int, op: OptionSmall) -> int:
	for NH in count(2):
		R = NH**2 * (2*NH - 1) / op.precision_ratio
		if not (R * M < 10**8 and R < 10**5 and L * R * M < 10**9):
			break
	return NH - 1

def create_family_vcf(family: Family, records: list[VCFFamilyRecord],
								is_mat_imputed: bool, num_families: int,
								ref_haps: list[list[int]],
								vcf: VCFSmall, op: OptionSmall) -> VCFImputable:
	NH = compute_upper_NH(family, len(vcf), num_families, op)
	if is_small(family, ref_haps, num_families, op):
		return VCFOneParentImputed(family.samples(), records, ref_haps,
											is_mat_imputed, op.map, vcf)
	elif is_small_ref(ref_haps, num_families, op):
		return VCFOneParentImputedRough(family.samples(), records, ref_haps,
													is_mat_imputed, op.map, vcf)
	elif NH >= 10:
		NH3 = min(20, NH)
		col = 1 if is_mat_imputed else 0
		gts = [ r.geno[col] for r in records ]
		filtered_ref_haps = filter_haplotypes(ref_haps, gts, NH3)
		return VCFOneParentImputedRough(family.samples(), records,
											filtered_ref_haps,
											is_mat_imputed, op.map, vcf)
	else:
		return VCFOneParentImputedFast(family.samples(), records,
											is_mat_imputed, op.map, vcf)

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
									ref_haps: list[list[int]],
									families: list[Family],
									non_imputed_parents: list[str],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFImputable] = []
	for family in families:
		vcf_family = VCFFamily.create_by_two_vcfs(imputed_vcf,
													orig_vcf, family.samples())
		is_mat_imputed = vcf_family.pat() in non_imputed_parents
		vcf = create_family_vcf(family, vcf_family.records, is_mat_imputed,
										len(families), ref_haps, orig_vcf, op)
		vcfs.append(vcf)
	
	if not vcfs:
		return None
	
	for vcf in vcfs:
		vcf.impute()
	
	print("%d families whose one parent is imputed and the other parent is"
									" known have been imputed." % len(vcfs))
	vcfs_base: list[VCFGenoBase] = [ vcf for vcf in vcfs ]
	new_vcf = VCFGeno.join(vcfs_base, orig_vcf.samples)
	return new_vcf

__all__ = ['impute', 'impute_family']
