from __future__ import annotations

# coding: utf-8
# SelfNonImputedFamily.py
# 自殖で一つもimputedが無い家系

from typing import Optional

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from Map import *
from KnownFamily import *
from VCFSelfImputable import VCFSelfImputable
from VCFSelfNoImputed import VCFSelfNoImputed
from VCFSelfNoImputedRough import VCFSelfNoImputedRough
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall


#################### SelfNonImputedFamily ####################

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

# is_smallが通るギリギリのNHを求める
def compute_upper_NH(family: Family, M: int, L: int, op: OptionSmall) -> int:
	for NH in count(2):
		R = NH**2 * (2*NH - 1) / op.precision_ratio
		if not (R * M < 10**8 and R < 10**5 and L * R * M < 10**9):
			break
	return NH - 1

def create_family_vcf(family: Family,
						records: list[GenoRecord],
						num_families: int,
						ref_haps: list[list[int]],
						vcf: VCFSmall,
						op: OptionSmall) -> Optional[VCFSelfImputable]:
	samples = [family.mat] + family.progenies
	NH = compute_upper_NH(family, len(vcf), num_families, op)
	if is_small(family, ref_haps, num_families, op):
		return VCFSelfNoImputed(samples, records, ref_haps, op.map, vcf)
	elif is_small_ref(ref_haps, num_families, op):
		return VCFSelfNoImputedRough(samples, records, ref_haps, op.map, vcf)
	elif NH >= 10:
		NH3 = min(20, NH)
		gts = [ r.geno[0] for r in records ]
		filtered_ref_haps = filter_haplotypes(ref_haps, gts, NH3)
		return VCFSelfNoImputedRough(samples, records,
										filtered_ref_haps, op.map, vcf)
	else:
		return None

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
										families: list[KnownFamily],
										op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		samples = [family.mat] + family.progenies
		vcf = VCFGeno.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
		vcf1 = create_family_vcf(family, vcf.records, len(families),
													ref_haps, orig_vcf, op)
		if vcf1 is not None:
			vcf1.impute()
			vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d self families have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
