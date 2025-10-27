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
from OptionSmall import OptionSmall

def is_small(ref_haps: list[list[int]], op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5

def impute(samples: list[str], orig_vcf: VCFSmall, ref_haps: list[list[int]],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	if not samples:
		return None
	
	vcf = VCFGeno.extract_samples(samples, orig_vcf)
	if is_small(ref_haps, op):
		vcf1 = VCFOrphan(samples, vcf.records, ref_haps, op.map, orig_vcf)
		vcf1.impute()
		print("%d orphan samples have been imputed." % len(samples))
		return vcf1
	else:
		return None

__all__ = ['impute']
