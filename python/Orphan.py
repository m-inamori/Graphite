from __future__ import annotations

# coding: utf-8
# Orphan.py
# 孤立したサンプルをimputeする

from typing import Optional

from VCFFamily import *
from Map import *
from VCFOrphan import VCFOrphan

def is_small(ref_haps: list[list[int]]) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1)
	return R * M < 10**8 and R < 10**5

def impute(samples: list[str], orig_vcf: VCFSmall,
			ref_haps: list[list[int]], gmap: Map) -> Optional[VCFSmallBase]:
	if not samples:
		return None
	
	vcf = orig_vcf.select_samples(samples)
	if is_small(ref_haps):
		vcf1 = VCFOrphan(vcf.header, vcf.records, ref_haps, gmap)
		vcf1.impute()
		print("%d orphan samples have been imputed." % len(samples))
		return vcf1
	else:
		return None

__all__ = ['impute']
