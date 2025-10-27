from __future__ import annotations

# coding: utf-8
# ProgenyImputedFamily.py
# 片親がknownで後代が一つ以上imputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFProgenyImputed import *
from Map import *
from KnownFamily import KnownFamily
from OptionSmall import OptionSmall

def impute(orig_vcf: VCFSmall, merged_vcf: VCFGenoBase,
			families: list[KnownFamily], imputed_progenies: list[list[str]],
			ref_haps: list[list[int]], op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for i in range(len(families)):
		family = families[i]
		parent = family.mat if family.mat_known else family.pat
		progeny = imputed_progenies[i][0]
		samples = [parent, progeny]
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf, orig_vcf, samples)
		vcf1 = VCFProgenyImputed(samples, vcf.records, ref_haps,
											family.mat_known, op.map, orig_vcf)
		vcf1.impute()
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose progeny is imputed have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
