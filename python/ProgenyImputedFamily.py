from __future__ import annotations

# coding: utf-8
# ProgenyImputedFamily.py
# 片親がknownで後代が一つ以上imputeされている家系を補完する

from VCFFamily import *
from VCFProgenyImputed import *
from Map import *
from KnownFamily import KnownFamily
from typing import Optional

def impute(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
			families: list[KnownFamily], imputed_progenies: list[list[str]],
			ref_haps: list[list[int]], gmap: Map) -> Optional[VCFSmall]:
	vcfs: list[VCFSmallBase] = []
	for i in range(len(families)):
		family = families[i]
		parent = family.mat if family.mat_known else family.pat
		progeny = imputed_progenies[i][0]
		samples = [parent, progeny]
		vcf = VCFSmall.create_by_two_vcfs(merged_vcf, orig_vcf, samples)
		vcf1 = VCFProgenyImputed(vcf.header, vcf.records, ref_haps,
													family.mat_known, gmap)
		vcf1.impute()
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose progeny is imputed have been imputed." %
															len(families))
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
