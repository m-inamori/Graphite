from __future__ import annotations

# coding: utf-8
# Orphan.py
# 孤立したサンプルをimputeする

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from Map import *
from VCFOrphan import VCFOrphan
from VCFOrphanRough import VCFOrphanRough
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall

def is_small(ref_haps: list[list[int]], op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5

# is_smallが通るギリギリのNHを求める
def compute_upper_NH(M: int, op: OptionSmall) -> int:
	for NH in count(2):
		R = NH**2 * (2*NH - 1) / op.precision_ratio
		if not (R * M < 10**8 and R < 10**5):
			break
	return NH - 1

def impute(samples: list[str], orig_vcf: VCFSmall, ref_haps: list[list[int]],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	if not samples:
		return None
	
	vcf = VCFGeno.extract_samples(samples, orig_vcf)
	NH = compute_upper_NH(len(vcf), op)
	if is_small(ref_haps, op):
		vcf1 = VCFOrphan(samples, vcf.records, ref_haps, op.map, orig_vcf)
		vcf1.impute()
		print("%d orphan samples have been imputed." % len(samples))
		return vcf1
	elif NH >= 10:
		NH2 = min(20, NH)
		ref_haps_table = []
		for i, sample in enumerate(samples):
			gts = [ r.geno[i] for r in vcf.records ]
			ref_haps_table.append(filter_haplotypes(ref_haps, gts, NH2))
		vcf2 = VCFOrphanRough(samples, vcf.records,
								ref_haps_table, op.map, orig_vcf)
		vcf2.impute()
		print("%d orphan samples have been imputed." % len(samples))
		return vcf2
	else:
		return None

__all__ = ['impute']
